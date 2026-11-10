"""Calibration packs: only fill what the input leaves open, and are traceable."""
import json

from sasnet import Network, JsonCalibration, Calibration
from sasnet.report import markdown
from sasnet import calibration as cal_mod


def _net(cd=None):
    lab = {"type": "labyrinth", "name": "seal", "a": "hp", "b": "lp",
           "clearance": 0.3e-3, "radius": 0.08, "teeth": 4, "pitch": 4e-3}
    if cd is not None:
        lab["cd"] = cd
    return {"nodes": [{"name": "hp", "boundary": True, "p": 4e5, "T": 600},
                      {"name": "lp", "boundary": True, "p": 1e5, "T": 600}],
            "elements": [lab]}


RULES = {"name": "test-pack", "version": "1", "rules":
         {"labyrinth": {"cd": {"c0": 0.80, "per": {"clearance_mm": -0.2}, "min": 0.5, "max": 0.9}}}}


def test_no_calibration_uses_defaults():
    net = Network.from_dict(_net())
    assert net.calibration is None
    assert net.elements[0].cd == 0.7


def test_json_rule_fills_missing_cd():
    net = Network.from_dict(_net(), RULES)
    assert abs(net.elements[0].cd - (0.80 - 0.2 * 0.3)) < 1e-12
    assert net.calibration.applied == [("seal", "cd", net.elements[0].cd)]


def test_explicit_value_wins():
    net = Network.from_dict(_net(cd=0.65), RULES)
    assert net.elements[0].cd == 0.65
    assert net.calibration.applied == []


def test_clip_to_bounds():
    rules = json.loads(json.dumps(RULES))
    rules["rules"]["labyrinth"]["cd"]["c0"] = 2.0
    net = Network.from_dict(_net(), rules)
    assert net.elements[0].cd == 0.9


def test_json_file_and_report(tmp_path):
    f = tmp_path / "pack.json"
    f.write_text(json.dumps(RULES))
    net = Network.from_dict(_net(), str(f))
    text = markdown(net, net.solve())
    assert "test-pack 1" in text and "Calibrated parameters" in text


def test_python_pack_subclass():
    class Pack(Calibration):
        def fill(self, kind, params):
            return {"cd": 0.61} if kind == "labyrinth" else {}
    net = Network.from_dict(_net(), Pack(name="py-pack"))
    assert net.elements[0].cd == 0.61


def test_unknown_pack_raises():
    try:
        cal_mod.load("no-such-pack")
    except LookupError:
        return
    assert False
