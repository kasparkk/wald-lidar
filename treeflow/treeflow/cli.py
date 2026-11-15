"""Command line:  python -m treeflow run INPUT --out RESULTS [--bbox x0 y0 x1 y1] [--ortho photo.tif]
                 [--reference trees_ref.csv] [--config settings.yaml]
INPUT can be one .las/.laz file or a folder of them."""
import argparse, os, sys
from .pipeline import run, run_folder, run_compare


def main(argv=None):
    p = argparse.ArgumentParser(prog='treeflow', description='Point cloud in, single trees out.')
    sub = p.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run', help='process a file or a folder')
    r.add_argument('input'); r.add_argument('--out', required=True)
    r.add_argument('--bbox', nargs=4, type=float, metavar=('X0', 'Y0', 'X1', 'Y1'))
    r.add_argument('--ortho', help='RGB+NIR aerial photo (GeoTIFF) in the same CRS')
    r.add_argument('--reference', help='CSV with x,y[,height_m,dbh_cm] of known trees')
    r.add_argument('--config', help='YAML file with settings (see config.example.yaml)')
    cp = sub.add_parser('compare', help='run several detectors on the same input')
    cp.add_argument('input'); cp.add_argument('--out', required=True)
    cp.add_argument('--detectors', nargs='+', default=['treeflow', 'lidr-dalponte', 'lidr-silva', 'lidr-li2012'])
    cp.add_argument('--bbox', nargs=4, type=float, metavar=('X0', 'Y0', 'X1', 'Y1'))
    cp.add_argument('--ortho'); cp.add_argument('--reference'); cp.add_argument('--config')
    v = sub.add_parser('visualise', help='already segmented trees -> 3D model bundle for the viewer')
    v.add_argument('input', help='point cloud with a tree id per point, or a folder of single-tree files')
    v.add_argument('--out', required=True); v.add_argument('--site'); v.add_argument('--title')
    v.add_argument('--tree-field', help='name of the tree id field (found automatically if common)')
    v.add_argument('--labels', help='CSV with column id plus values to show (species, source, ...)')
    v.add_argument('--tile', type=float, default=100.0); v.add_argument('--config')
    v.add_argument('--ortho'); v.add_argument('--reference')
    a = p.parse_args(argv)
    cfg = {}
    if a.config:
        import yaml
        cfg = yaml.safe_load(open(a.config)) or {}
    kw = dict(cfg=cfg, ortho=a.ortho, reference=a.reference)
    if a.cmd == 'visualise':
        from .visualise import visualise
        info = visualise(a.input, a.out, site=a.site, title=a.title, tree_field=a.tree_field,
                         labels=a.labels, tile=a.tile, cfg=cfg)
        print(f"{info['trees']} tree models in {len(info['tiles'])} tiles ({info['trees_from']}) -> {a.out}/index.json")
        return
    if a.cmd == 'compare':
        t = run_compare(a.input, a.out, a.detectors, bbox=a.bbox, **kw)
        for r in t:
            print({k: r[k] for k in r if k in ('detector', 'trees', 'trees_per_ha', 'found_pct', 'correct_pct', 'height_rmse_m')})
        return
    if os.path.isdir(a.input):
        ov = run_folder(a.input, a.out, **kw)
        print(f'{len(ov)} files processed → {a.out}/overview.csv')
    else:
        rows, s, v = run(a.input, a.out, bbox=a.bbox, **kw)
        print(f"{s['trees']} trees ({s.get('trees_per_ha')} per ha), method: {s.get('method')} → {a.out}/report.html")
        if v:
            print('check:', v)


if __name__ == '__main__':
    sys.exit(main())
