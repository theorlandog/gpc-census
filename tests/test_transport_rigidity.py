import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.transport_rigidity import audit  # noqa: E402

def test_observable_5_2_product_is_hopping_rigid():
    r=audit(5,2,5,2)
    assert r['pair_incidence_full_column_rank']
    assert r['carrier_preserving_hopping_dimension']==0

def test_4_2_product_has_hopping_but_loses_pair_rank():
    r=audit(4,2,4,2)
    assert r['carrier_preserving_hopping_dimension']==16
    assert r['pair_incidence_rank']==11
    assert not r['pair_incidence_full_column_rank']
