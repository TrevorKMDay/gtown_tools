import argparse as ap
import tempfile as tf
import subprocess as sp
import os

from pprint import pp

from PIL import Image
from nilearn.image import resample_to_img

# Arguments =====

parser = ap.ArgumentParser()

parser.add_argument("--xyz", required=True, nargs=3, metavar=("X", "Y", "Z"),
                    help="x/y/z coordinates, in voxel space.")

parser.add_argument("-t", "--threshold", nargs=2,
                    help="Threshold: lower, upper")

parser.add_argument("-o", "--output_dir",
                    help="Directory to save files to.")

parser.add_argument("-L", "--label", action="store_true",
                    help="Include slice label (option: -L)")

parser.add_argument("-H", "--height", metavar="inches", default=1.5,
                    help="Height to make images, in inches (300 dpi); "
                         "default=1.5")


parser.add_argument("stat_image", nargs="+",
                    help="Statistic image(s)")

args = parser.parse_args()

pp(args)

# Create image =====

xyz = ["x", "y", "z"]

def create_underlay(stat_image):

    print("Resampling T1 to requested space")

    fsl_mni = "/Users/tkmd/bin/fsl/data/standard/MNI152_T1_0.5mm.nii.gz"
    t1_underlay = resample_to_img(fsl_mni, stat_image)

    return(t1_underlay)

def create_image(coords, threshold, stat_image, new_height = 2, label = False, 
                 o_dir=None):

    if o_dir is None:
        o_dir = os.path.dirname(stat_image)

    stat_image_bn = os.path.basename(stat_image).rsplit(".")[0]
    output_file = os.path.join(o_dir, f"{stat_image_bn}.png")

    with tf.TemporaryDirectory(delete=False) as td:

        files = [os.path.join(td, f"{x}.png") for x in xyz]

        # Create overlay

        underlay = create_underlay(stat_image)
        underlay_f = os.path.join(td, "t1.nii.gz")
        underlay.to_filename(underlay_f)

        overlay_f = os.path.join(td, "overlay.nii.gz")
        overlay_cmd = ["overlay", "0", "0", 
                       underlay_f, "-a", 
                       stat_image, threshold[0], threshold[1],
                       overlay_f]
        
        # print(overlay_cmd)
        
        sp.run(overlay_cmd)

        label_str = "-L" if label else ""

        # Run slicer
        cmd = ["slicer", 
                    label_str,
                    overlay_f,
                    "-l", "Render1",
                    f"-x", f"-{coords[0]}", files[0], 
                    f"-y", f"-{coords[1]}", files[1], 
                    f"-z", f"-{coords[2]}", files[2], 
                    "-c"
                ]

        # print(cmd)
        sp.run(cmd)
            
        ims = [Image.open(x) for x in files]
        
        xyz_image = Image.new('RGB', 
                              (sum([ims[0].width, ims[1].width, ims[2].width]),
                               max([ims[0].height, ims[1].height, 
                                    ims[2].height])
                              ))

        # Concatenate images horizontally
        xyz_image.paste(ims[0], (0, 0))
        xyz_image.paste(ims[1], (ims[0].width, 0))
        xyz_image.paste(ims[2], (ims[0].width + ims[1].width, 0))

        # Resize image to 600 px tall (arbitrary)
        height_px = new_height * 300
        ratio = height_px / xyz_image.size[0] 
        new_size = (int(float(xyz_image.size[0]) * ratio),
                    int(float(xyz_image.size[1]) * ratio))
        
        xyz_image = xyz_image.resize(new_size, resample=Image.NEAREST)

        xyz_image.save(output_file)
        print(output_file)

# Main loop =====

new_height = float(args.height)

for stat_image in args.stat_image:

    create_image(args.xyz, args.threshold, stat_image, new_height, 
                 args.label, args.output_dir)
