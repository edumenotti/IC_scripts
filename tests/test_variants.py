import pytest

from mutation_hunter import variants

SEQS = {"S1": "MKA*", "S2": "MKA*", "S3": "MKA*", "S4": "MKV*", "S5": "MXA*"}


def test_parse_query():
    assert variants.parse_query("YNL079C") == ("YNL079C", None)
    assert variants.parse_query("YAR019W-A:A123*") == ("YAR019W-A", ("A", 123, "*"))
    with pytest.raises(ValueError, match="bad query"):
        variants.parse_query("YNL079C:123V")


def test_rare_variants_uses_major_allele_and_ignores_X():
    rv = variants.rare_variants("G", SEQS, max_freq=0.2, min_carriers=1)
    assert rv["variant"].tolist() == ["A3V"]
    assert rv.loc[0, "carriers"] == ["S4"]
    assert variants.rare_variants("G", SEQS, max_freq=0.1, min_carriers=1).empty  # 1/5 = 20% > 10%


def test_carriers_of_validates():
    assert variants.carriers_of("G", SEQS, "A", 3, "V") == ["S4"]
    with pytest.raises(ValueError, match="major allele is A, not V"):
        variants.carriers_of("G", SEQS, "V", 3, "A")
    with pytest.raises(ValueError, match="out of range"):
        variants.carriers_of("G", SEQS, "A", 99, "V")
    with pytest.raises(ValueError, match="no strain carries"):
        variants.carriers_of("G", SEQS, "A", 3, "W")


def test_carrier_matrix_rows_are_strains_with_the_gene():
    G = variants.carrier_matrix(SEQS, {"A3V": ["S4"]})
    assert G["A3V"].to_dict() == {"S1": 0.0, "S2": 0.0, "S3": 0.0, "S4": 1.0, "S5": 0.0}


def test_unaligned_input_is_rejected():
    with pytest.raises(ValueError, match="not aligned"):
        variants.rare_variants("G", {"S1": "MK", "S2": "MKA"}, 0.5, 1)
