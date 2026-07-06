import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import calibrate_genome_track_margin as cm

def test_locus_join():
    h = ">lcl|NC_088360.1_cds_XP_078753709.1_1 [gene=L] [location=join(4611..4721,5584..5704,6641..6776)] [gbkey=CDS]"
    assert cm.parse_cds_true_locus(h) == ("NC_088360.1", 4611, 6776)

def test_locus_complement_partial():
    h = ">lcl|NC_088360.1_cds_XP_1.1_2 [location=complement(join(40140..40160,48651..>48764))]"
    assert cm.parse_cds_true_locus(h) == ("NC_088360.1", 40140, 48764)

def test_locus_single_exon():
    h = ">lcl|NW_999.1_cds_XP_9.1_3 [location=100..250]"
    assert cm.parse_cds_true_locus(h) == ("NW_999.1", 100, 250)

def test_locus_missing_fields():
    with pytest.raises(ValueError):
        cm.parse_cds_true_locus(">lcl|NC_1.1_cds_XP_9.1_1 [gene=x]")   # no [location=]
    with pytest.raises(ValueError):
        cm.parse_cds_true_locus(">bad_header_no_cds [location=1..2]")   # no _cds_
