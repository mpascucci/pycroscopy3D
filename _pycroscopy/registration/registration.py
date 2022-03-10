import ants 
import multipagetiff as mtif

def load_for_ants(img_path):
    """Load an ANTs image with the right axis order"""
    stack = mtif.read_stack(img_path)
    return ants.from_numpy(stack.pages.astype(float))

def register_with_ANTs(paths, template_path, out_folder, type_of_transform="Rigid"):
    template = load_for_ants(template_path)

    for path in tqdm(paths):
        out_path = os.path.join(out_folder, os.path.basename(path))
        to_register = load_for_ants(path)
    
        areg = ants.registration(fixed=template, moving=to_register, 
                                  type_of_transform=type_of_transform)
        corrected = areg['warpedmovout'].numpy()
        stack = mtif.Stack(corrected)
        stack._dtype_out = np.uint16
        stack._apply_normalization()
        mtif.write_stack(stack, out_path)