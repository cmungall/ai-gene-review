import pytest

@pytest.mark.parametrize("a, b, c", [(1, 2, 3), (4, 5, 9)])
def test_simple(a, b, c):
    assert a + b == c

