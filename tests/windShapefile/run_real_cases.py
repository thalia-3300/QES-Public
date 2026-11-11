#!/usr/bin/env python3
"""Reproducible Testing/Clare smoke tests; private GIS inputs are not vendored.

Uses only Python's standard library. All generated XML, logs and NetCDF outputs
are kept in a fresh temporary directory printed at startup. Height, roughness
and stability here are test assumptions, not recovered WindNinja metadata.
"""
import argparse
import copy
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def fields(parent, values):
    for name, value in values.items():
        ET.SubElement(parent, name).text = str(value)


def case(shp, clare=False, dem=None, building=False):
    root = ET.Element("QESWindsParameters")
    sim = ET.SubElement(root, "simulationParameters")
    fields(sim, dict(domain="48 48 30" if clare else "44 44 30",
                     cellSize="20 20 10", totalTimeIncrements=1,
                     maxIterations=1000, tolerance="1e-7", meshTypeFlag=0))
    if dem:
        fields(sim, dict(DEM=dem, originFlag=0, DEMDistancex=20, DEMDistancey=20))
    met = ET.SubElement(root, "metParams")
    fields(met, dict(z0_domain_flag=0))
    wind = ET.SubElement(met, "windShapefile")
    wind.append(ET.Comment("TEST ASSUMPTIONS: 10 m AGL, z0=0.1 m, neutral QES profile"))
    fields(wind, dict(file=shp, speedField="speed", directionField="dir",
                      speedUnits="m/s", directionUnits="degrees",
                      originX=277020 if clare else 374610,
                      originY=6252020 if clare else 5835000,
                      height=10, siteZ0=0.1, reciprocal=0, boundaryLayerFlag=1,
                      timeStamp="2026-09-17T16:30:00" if clare else "2026-09-10T21:09:00"))
    fields(ET.SubElement(root, "fileOptions"), dict(outputFlag=1, outputFields="all"))
    fields(ET.SubElement(root, "turbParams"), dict(method=0))
    if building:
        buildings = ET.SubElement(root, "buildingsParams")
        fields(buildings, dict(wallRoughness=0.001, rooftopFlag=1, upwindCavityFlag=2,
                              streetCanyonFlag=1, streetIntersectionFlag=0, wakeFlag=2,
                              highRiseFlag=0, sidewallFlag=1))
        fields(ET.SubElement(buildings, "rectangularBuilding"),
               dict(height=40, baseHeight=0, xStart=360, yStart=400,
                    length=80, width=80, buildingRotation=0))
    return root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--testing", type=Path, required=True)
    parser.add_argument("--clare", type=Path, required=True)
    parser.add_argument("--clare-dem", type=Path, required=True)
    parser.add_argument("--solver", choices=[1, 3, 4], type=int, default=1)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--require-convergence", action="store_true")
    args = parser.parse_args()
    build = args.build.resolve()
    binary = build / "qesWinds/qesWinds"
    checker = build / "tests/windShapefile/windShapefileTests"
    work = Path(tempfile.mkdtemp(prefix=f"qes-windninja-s{args.solver}-"))
    print(f"Results: {work}", flush=True)
    env = dict(os.environ, OMP_NUM_THREADS="4")

    def run(name, root, error=None):
        if not error:
            root.find("simulationParameters/maxIterations").text=str(args.iterations)
        xml = work / f"{name}.xml"
        ET.indent(root)
        ET.ElementTree(root).write(xml, encoding="utf-8", xml_declaration=True)
        command = [str(binary), "-q", str(xml), "-s", str(args.solver)]
        if not error:
            command += ["-t", "-w", "-o", str(work / name)]
        result = subprocess.run(command, cwd=work, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        (work / f"{name}.log").write_text(result.stdout)
        if error:
            if result.returncode == 0 or error not in result.stdout:
                raise RuntimeError(f"{name}: expected rejection '{error}'; see log")
        else:
            if result.returncode:
                raise RuntimeError(f"{name}: exit {result.returncode}; see log")
            if args.require_convergence:
                if "[Solver] CONVERGED" not in result.stdout or "NOT CONVERGED" in result.stdout:
                    raise RuntimeError(f"{name}: convergence criterion not met; see log")
                subprocess.run([str(checker), "--check-mass", str(work/f"{name}_windsWk.nc")],
                               check=True, env=env)
            subprocess.run([str(checker), "--check-output", str(work/f"{name}_windsOut.nc")],
                           check=True, env=env)
            subprocess.run([str(checker), "--check-turb", str(work/f"{name}_turbOut.nc")],
                           check=True, env=env)
        print(f"PASS {name}", flush=True)

    for shp, n in [(args.testing,142),(args.clare,51)]:
        subprocess.run([str(checker),"--real",str(shp.resolve()),str(n)],check=True,env=env)
    testing = case(str(args.testing.resolve()))
    run("testing_flat", testing)
    run("testing_building", case(str(args.testing.resolve()),building=True))
    run("clare_flat", case(str(args.clare.resolve()),clare=True))
    run("clare_terrain", case(str(args.clare.resolve()),clare=True,dem=str(args.clare_dem.resolve())))
    # Keep the original sensor path working when windShapefile is absent.
    legacy = copy.deepcopy(testing)
    met = legacy.find("metParams")
    met.remove(met.find("windShapefile"))
    sensor = ET.SubElement(met,"sensor")
    fields(sensor,dict(site_coord_flag=1,site_xcoord=440,site_ycoord=440))
    fields(ET.SubElement(sensor,"timeSeries"),
           dict(timeStamp="2026-09-10T21:09:00",boundaryLayerFlag=1,
                siteZ0=.1,reciprocal=0,height=10,speed=5,direction=270))
    run("legacy_sensor",legacy)

    # End-to-end XML validation, in addition to the synthetic C++ grid tests.
    for name, path, value, error in [
        ("uncovered", "metParams/windShapefile/originX", "374000", "coverage"),
        ("multiple_times", "simulationParameters/totalTimeIncrements", "2", "one timestep"),
        ("unknown_units", "metParams/windShapefile/speedUnits", "unknown", "speedUnits"),
        ("missing_field", "metParams/windShapefile/speedField", "absent", "attribute field"),
    ]:
        invalid = copy.deepcopy(testing)
        invalid.find(path).text = value
        run(name,invalid,error)
    for name, tag, value, error in [
        ("halo", "halo_x", 20, "zero halo"),
        ("rotated", "domainRotation", 10, "zero domainRotation"),
    ]:
        invalid = copy.deepcopy(testing)
        fields(invalid.find("simulationParameters"),{tag:value})
        run(name,invalid,error)
    invalid = copy.deepcopy(testing)
    wind = invalid.find("metParams/windShapefile")
    wind.remove(wind.find("height"))
    run("missing_height",invalid,"height")
    invalid = copy.deepcopy(testing)
    fields(invalid.find("metParams/windShapefile"),{"timesteps":2})
    run("unknown_option",invalid,"unknown or repeated option")
    invalid = copy.deepcopy(testing)
    fields(invalid.find("metParams/windShapefile"),{"height":20})
    run("repeated_height",invalid,"unknown or repeated option")
    invalid = copy.deepcopy(testing)
    met = invalid.find("metParams")
    met.append(copy.deepcopy(met.find("windShapefile")))
    run("repeated_source",invalid,"Only one windShapefile")
    invalid = copy.deepcopy(testing)
    fields(invalid.find("metParams"),{"sensorName":"unused.xml"})
    run("mixed_sources",invalid,"cannot be combined")
    invalid = case(str(args.clare.resolve()),clare=True,dem=str(args.clare_dem.resolve()))
    invalid.find("metParams/windShapefile/originX").text="277021"
    run("misaligned_dem",invalid,"DEM crop origin")
    print(f"All real-case smoke tests passed. Logs/XML/outputs: {work}",flush=True)


if __name__ == "__main__":
    main()
