from trebouxia_pipeline.cli import parse_ann_field


def test_parse_ann_field_basic():
    info = "DP=20;ANN=T|missense_variant|MODERATE|geneA|g0001|transcript|tx1|protein_coding|1/1|c.10A>T|p.K4M|10/100|10/90|4/30||"
    records = parse_ann_field(info)
    assert len(records) == 1
    assert records[0]["effect"] == "missense_variant"
    assert records[0]["impact"] == "MODERATE"
    assert records[0]["gene_id"] == "g0001"
