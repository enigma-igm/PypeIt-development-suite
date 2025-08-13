from astroquery.mast import MastMissionsClass

# --- setup (use your real token locally, not in shared code) ---
token = "9465311c9b2544afaa8e8c1917f63f5a"
MastClass = MastMissionsClass(mission='JWST')
MastClass.login(token=token)

program = '9180'  # '9180' = fixed-slit (DW); '4713' = MSA (Masquerade)
#program = '4713'  # '9180' = fixed-slit (DW); '4713' = MSA (Masquerade)


# --- program-specific config (keeps your originals) ---
if program == '9180':
    exptype_sci = 'NRS_FIXEDSLIT'     # FS science
    download_dir = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/9180'
    do_ta = True                      # only FS gets TA
elif program == '4713':
    exptype_sci = 'NRS_MSASPEC'       # MSA science
    download_dir = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/4713'
    do_ta = False
else:
    raise ValueError(f"Unknown program: {program}")

# ----------------------------
# 0) Full program inventory
# ----------------------------
all_rows = MastClass.query_criteria(program=program)

# ----------------------------
# 1) SCIENCE (like before)
# ----------------------------
datasets_sci = MastClass.query_criteria(
    program=program,
    exp_type=exptype_sci,
    productLevel='1b'  # lean: Level 1b only
)

products_sci = MastClass.get_unique_product_list(datasets_sci) if len(datasets_sci) > 0 else []

if exptype_sci == 'NRS_FIXEDSLIT':
    filtered_sci = MastClass.filter_products(
        products_sci, file_suffix=['_uncal'], extension='fits'
    )
else:
    filtered_sci = MastClass.filter_products(
        products_sci, file_suffix=['_uncal', '_rate', '_msa'], extension='fits'
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
