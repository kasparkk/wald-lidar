"""DetailView (Schindler & Frey, Uni Freiburg) embeddings for single-tree files: 256-d hidden layer + 33 class probabilities,
averaged over augmentations. Needs a clone of github.com/JulFrey/DetailView (patched: densenet weights=None,
map_location=cpu, read_las in metres) and the weights from FreiData (doi 10.60493/f850a-bb152).
  DETAILVIEW_DIR=... DETAILVIEW_WEIGHTS=... python3 dv_extract.py trees.csv LAS_FOLDER/ out.csv 3
trees.csv: filename,species_id,tree_H (species_id -999 if unknown)."""
import sys, os, numpy as np, pandas as pd, torch
DV = os.environ.get('DETAILVIEW_DIR', 'DetailView'); sys.path.insert(0, DV); os.chdir(DV)
torch.set_grad_enabled(False); torch.set_num_threads(2)
import parallel_densenet as net
from torchvision import transforms
csv, folder, out, naug = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
model = net.SimpleView(n_classes=33, n_views=7)
model.load_state_dict(torch.load(os.environ.get('DETAILVIEW_WEIGHTS', 'model_ft_202412171652_3'), map_location='cpu')); model.eval()
H = {}
model.classifier[1].register_forward_hook(lambda m, i, o: H.__setitem__('z', o.detach().numpy().copy()))
ds = net.TrainDataset_AllChannels(csv, folder, img_trans=transforms.Compose([transforms.RandomVerticalFlip(0.5)]), pc_rotate=True,
    height_noise=0.01, test=True, res=256, n_sides=4, height_mean=15.2046, height_sd=9.5494, projection_backend='torch')
dl = torch.utils.data.DataLoader(ds, batch_size=1, shuffle=False, num_workers=0)
Z, P, names = {}, {}, []
for a in range(naug):
    for x, h, path in dl:
        p = torch.softmax(model(x, h), 1).numpy()[0]; n = path[0]
        Z.setdefault(n, []).append(H['z'][0]); P.setdefault(n, []).append(p)
    print('aug', a + 1, flush=True)
rows = [dict(file=n, **{f'z{i}': v for i, v in enumerate(np.mean(Z[n], 0))}, **{f'p{i}': v for i, v in enumerate(np.mean(P[n], 0))}) for n in Z]
pd.DataFrame(rows).to_csv(out, index=False)
