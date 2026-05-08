import os
from astroquery.mast import MastMissionsClass


# --- setup (use your real token locally, not in shared code) ---
token = "e5d75cc5337e4dbaa0f6784fdf4480bc" # Jiamu's token on Jan12, 2026
MastClass = MastMissionsClass(mission='JWST')
MastClass.login(token=token)

program = '5664' # JWST IFU GO program for LRDs


# --- program-specific config (keeps your originals) ---
if program == '9180' or program == '3417' or program == '7491':
    exptype_sci = 'NRS_FIXEDSLIT'     # FS science
    download_dir = os.path.join('/Users/joe/jwst_redux/Raw/NIRSPEC_FS/', program)
    do_ta = True                      # only FS gets TA
elif program == '4713':
    exptype_sci = 'NRS_MSASPEC'       # MSA science
    download_dir = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/4713'
    do_ta = False
elif program == '5664':
    exptype_sci = 'NRS_IFU'       # IFU science
    download_dir = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_IFU/5664/'
    do_ta = False
else:
    raise ValueError(f"Unknown program: {program}")

# ----------------------------
# 0) Full program inventory
# ----------------------------
all_rows = MastClass.query_criteria(program=program)

datasets_sci = all_rows 
#datasets_sci = MastClass.query_criteria(
#    program=program,
#    exp_type=exptype_sci)
# Removing produclevel restriction to download things that may have problems
#    productLevel='1b'  # lean: Level 1b only
#)

products_sci = MastClass.get_unique_product_list(datasets_sci) if len(datasets_sci) > 0 else []

if exptype_sci == 'NRS_FIXEDSLIT':
    filtered_sci = MastClass.filter_products(
        products_sci, file_suffix=['_uncal', '_rate'], extension='fits'
    )
else:
    # restrict products to one pointing
    # target_obs = 'jw05664002001_02103' # For GS-13971 which is the only one available now, also jw05664002001_02101
    target_obs = 'jw05664006001_04101' # Try GN-15498
    mask = [target_obs in fn for fn in products_sci['filename']]
    products_sci_subset = products_sci[mask]
    filtered_sci = MastClass.filter_products(
        products_sci_subset, file_suffix=['_cal'], extension='fits'
    )

if len(filtered_sci) > 0:
    MastClass.download_products(filtered_sci, download_dir=download_dir, verbose=True)
else:
    print(f"[SCI] No products matched for program {program}.")

# ---------------------------------------------
# 2) TA for FIXED-SLIT only (client-side filter)
# ---------------------------------------------
if do_ta:
    # A) Primary: filter locally by exp_type == 'NRS_WATA'
    # (works even when the server-side exp_type filter returns nothing)
    try:
        ta_rows = all_rows[all_rows['exp_type'] == 'NRS_WATA']
    except Exception:
        ta_rows = []


    if len(ta_rows) == 0:
        print(f"[TA] No TA rows found in program {program} (no NRS_WATA or intentType=acquisition).")
    else:
        # IMPORTANT: pass the filtered table directly; don't re-query with exp_type
        prods_wata = MastClass.get_unique_product_list(ta_rows)

        keep_wata = MastClass.filter_products(
            prods_wata,
            file_suffix=['_uncal', '_rate', '_cal', '_i2d'],  # raw + processed for WATA
            extension='fits'
        )

        if len(keep_wata) > 0:
            MastClass.download_products(keep_wata, download_dir=download_dir, verbose=True)
        else:
            print(f"[TA] TA rows exist, but no matching WATA products with desired suffixes in program {program}.")
