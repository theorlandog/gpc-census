%global srcname gpc_census

Name:           gpc-census
# `make rpm` rewrites this line with the version from pyproject.toml
Version:        0.1.0
Release:        1%{?dist}
Summary:        Exact extremal states for fermionic natural-occupation-number polytopes

License:        GPL-3.0-or-later
URL:            https://github.com/theorlandog/gpc-census
Source0:        %{srcname}-%{version}.tar.gz

BuildArch:      noarch
BuildRequires:  python3-devel
BuildRequires:  pyproject-rpm-macros

%global _description %{expand:
gpc-census implements an algorithmic pipeline that constructs exact
extremal states for fermionic natural-occupation-number (moment)
polytopes. It can be used as a Python library (import gpc_census) or
through the gpc-census command-line tool.}

Recommends: python3dist(ortools)

%description %_description

%prep
%autosetup -n %{srcname}-%{version}

%generate_buildrequires
%pyproject_buildrequires -g test

%build
%pyproject_wheel

%install
%pyproject_install
%pyproject_save_files -l %{srcname}

%check
%pytest

%files -f %{pyproject_files}
%doc README.md
%{_bindir}/gpc-census

%changelog
* Fri Jul 17 2026 Jamie Orlando <jamie@orlandonh.com> - 0.1.0-1
- Initial package
