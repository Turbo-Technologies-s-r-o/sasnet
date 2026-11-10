"""Command line:

    python -m sasnet network.json [report.md] [--calibration PACK]
    python -m sasnet --list-calibrations

PACK is the name of an installed calibration pack or a .json rules file.
Without it the first-principles defaults are used (or the "calibration"
field of the network file, if present).
"""
import argparse
import sys

from .network import Network
from .report import markdown
from . import calibration


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m sasnet", description="SAS network solver")
    ap.add_argument("network", nargs="?")
    ap.add_argument("report", nargs="?")
    ap.add_argument("--calibration", "-c", default=None)
    ap.add_argument("--list-calibrations", action="store_true")
    a = ap.parse_args(argv)
    if a.list_calibrations:
        print("\n".join(calibration.available()) or "no calibration packs installed")
        return 0
    if not a.network:
        print(__doc__)
        return 2
    net = Network.from_json(a.network, a.calibration)
    r = net.solve()
    text = markdown(net, r, title=a.network)
    if a.report:
        with open(a.report, "w", encoding="utf-8") as fh:
            fh.write(text)
    print(text)
    return 0 if r.converged else 1


if __name__ == "__main__":
    sys.exit(main())
