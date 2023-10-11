import numpy as np
from ldpc import bposd_decoder
from bposd.hgp import hgp
import panqec.codes
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="distance", type=int)
parser.add_argument("-p", help="error rate times 10000", type=int)
parser.add_argument("-o", help="output file to write to", type=str)
args = parser.parse_args()

d = args.d
p = args.p / 10000.
c = panqec.codes.surface_3d.Planar3DCode(d,d,d)
H = c.Hx
L = c.logicals_x[0][0:c.n]

if d == 3:
    order = 19
else:
    order = 60
bpd=bposd_decoder(
    H,
    error_rate=p,
    max_iter=c.n, 
    bp_method="ms",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    osd_order=order #the osd search depth
)


succ = 0
fail = 0
for _ in range(1000000):
    error = (np.random.uniform(size=(c.n)) < p).astype(int)
    syndrome=H@error %2
    bpd.decode(syndrome)

    #Decoding is successful if the residual error commutes with the logical operators
    residual_error=(bpd.osdw_decoding+error) %2
    a=((L@residual_error)%2).any()
    if a: fail += 1
    else: succ += 1
    if _ % 100 == 0:
        print(_, fail / (fail+succ))
        # write result to disk
        if args.o is not None:
            with open(args.o, 'w') as fd:
                fd.write("{},{}\n".format(succ,fail))
