"""Spot-resilient checkpointing shared by the campaign scripts.

Each long campaign appends its results to a local output file and, on a fixed
wall-clock interval, snapshots that file to durable storage so a spot-instance
kill loses at most one interval of work. On a fresh instance the same file is
restored from durable storage before the run resumes, so a redeploy picks up
exactly where the previous instance stopped (the campaigns already skip
already-recorded work once the file is back in place).

Durable storage is the local disk by default: the campaigns write and flush
their output as they go, and each checkpoint additionally fsyncs it so a
same-disk restart is consistent. Passing an S3 bucket also mirrors the file to
S3, which is what actually survives losing the instance's ephemeral disk.
Credentials come from the named profile when one is given, otherwise from the
default boto3 credential chain (instance role, environment, shared config...).

boto3 is imported lazily, so disk-only runs need nothing beyond the standard
library. Install the extra with `uv pip install .[s3]` (or `boto3`) to enable
the S3 mirror.

Typical wiring:

    ck = Checkpointer.from_opts(default_out, opts)   # None if no output file
    if ck:
        ck.restore()                                 # pull prior progress
    with open(out_path, "a") as fh:
        ck and ck.attach(fh)
        ck and ck.install_signal_handlers()          # last-gasp on SIGTERM
        for task in todo:
            emit(fh, work(task))
            ck and ck.checkpoint()                    # honors the interval
        ck and ck.checkpoint(force=True)              # final flush + upload
"""
from __future__ import annotations

import os
import pathlib
import signal
import sys
import time

# Flags every checkpoint-capable script accepts; each consumes one value.
FLAGS = ("--out", "--checkpoint-interval", "--s3-bucket", "--s3-key", "--s3-profile")

DEFAULT_INTERVAL = 300.0


def extract_opts(argv):
    """Split checkpoint flags out of a positional argv.

    Returns (opts, remaining) where opts maps the de-dashed flag name to its
    string value and remaining is argv with those flag/value pairs removed, so
    a script's existing positional parsing keeps working untouched.
    """
    opts = {}
    remaining = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in FLAGS and i + 1 < len(argv):
            opts[a[2:].replace("-", "_")] = argv[i + 1]
            i += 2
        else:
            remaining.append(a)
            i += 1
    return opts, remaining


def add_arguments(parser, default_out=None):
    """Register the checkpoint flags on an argparse parser."""
    parser.add_argument("--out", default=str(default_out) if default_out else None,
                        help="output file to checkpoint (default: %(default)s)")
    parser.add_argument("--checkpoint-interval", type=float, default=DEFAULT_INTERVAL,
                        help="seconds between checkpoints (default: %(default)s)")
    parser.add_argument("--s3-bucket", default=None,
                        help="mirror the checkpoint to this S3 bucket")
    parser.add_argument("--s3-key", default=None,
                        help="S3 object key (default: the output file name)")
    parser.add_argument("--s3-profile", default=None,
                        help="AWS profile; omit to use the default credential chain")


def _s3_client(profile):
    import boto3  # lazy: only needed when a bucket is configured

    session = boto3.Session(profile_name=profile) if profile else boto3.Session()
    return session.client("s3")


class Checkpointer:
    def __init__(self, path, *, interval=DEFAULT_INTERVAL, bucket=None, key=None,
                 profile=None):
        self.path = pathlib.Path(path)
        self.interval = float(interval)
        self.bucket = bucket
        self.key = key or self.path.name
        self.profile = profile
        self._fh = None
        self._deadline = time.monotonic() + self.interval
        self._client = _s3_client(profile) if bucket else None

    @classmethod
    def from_opts(cls, default_path, opts):
        """Build from an extract_opts() dict; None when there is no output file."""
        path = opts.get("out", default_path)
        if not path:
            return None
        interval = opts.get("checkpoint_interval", DEFAULT_INTERVAL)
        return cls(path, interval=float(interval), bucket=opts.get("s3_bucket"),
                   key=opts.get("s3_key"), profile=opts.get("s3_profile"))

    def attach(self, fh):
        """Remember an open file handle so signal-driven checkpoints can fsync it."""
        self._fh = fh

    def restore(self):
        """If the output file is absent locally, pull the last S3 checkpoint.

        Returns True when a checkpoint was restored. A missing S3 object (a
        brand-new campaign) is not an error.
        """
        if self.path.exists() or not self._client:
            return False
        from botocore.exceptions import ClientError

        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self._client.download_file(self.bucket, self.key, str(self.path))
        except ClientError as e:
            if e.response.get("Error", {}).get("Code") in ("404", "NoSuchKey"):
                return False
            raise
        size = self.path.stat().st_size
        print(f"checkpoint: restored s3://{self.bucket}/{self.key} "
              f"-> {self.path} ({size} bytes)", file=sys.stderr, flush=True)
        return True

    def checkpoint(self, *, force=False):
        """Fsync the output and mirror to S3 when the interval has elapsed.

        Cheap to call in a tight loop: it returns immediately until the
        interval passes (or force=True), so callers need no timing logic.
        """
        now = time.monotonic()
        if not force and now < self._deadline:
            return False
        self._deadline = now + self.interval
        self._fsync()
        self._upload()
        return True

    def _fsync(self):
        if self._fh is None:
            return
        self._fh.flush()
        os.fsync(self._fh.fileno())

    def _upload(self):
        if not self._client or not self.path.exists():
            return
        # S3 object writes are atomic (a multipart upload only becomes visible
        # once completed), so uploading straight to the key is safe even while
        # the campaign keeps appending: a checkpoint may trail the newest lines
        # but is never a torn object.
        self._client.upload_file(str(self.path), self.bucket, self.key)
        print(f"checkpoint: uploaded {self.path} -> s3://{self.bucket}/{self.key}",
              file=sys.stderr, flush=True)

    def install_signal_handlers(self):
        """Force a final checkpoint on SIGTERM/SIGINT, then exit.

        Spot interruptions arrive as SIGTERM (typically ~2 min before the
        instance is reclaimed); flushing here saves the tail of the interval.
        """
        def handler(signum, _frame):
            print(f"checkpoint: caught signal {signum}, flushing final checkpoint",
                  file=sys.stderr, flush=True)
            try:
                self.checkpoint(force=True)
            finally:
                sys.exit(128 + signum)

        for name in ("SIGTERM", "SIGINT"):
            sig = getattr(signal, name, None)
            if sig is not None:
                try:
                    signal.signal(sig, handler)
                except (ValueError, OSError):
                    pass  # not in the main thread, or unsupported on this platform
