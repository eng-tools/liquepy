import pandas as pd
import glob
import os
import inspect


D_STR = "Depth (m)"
QC_STR = "qc (MPa)"
FS_STR = "fs (MPa)"
U2_STR = "u (MPa)"


def convert_file(ffp, out_fp, verbose=0):
    """
    Converts a CPT from the NZGD into the standard liquepy format.

    This solves the problem of dealing with multiple different file formats by converting to a single format.

    The additional metadata is mostly maintained in the file.

    Algorithm explained
    -------------------

    Algorithm works by trial and error:
     - Many different converters have been written to convert from an existing format to the liquepy format
     - This algorithm cycles through each converter and tries to recognise the existing file format
        and then apply the conversion
     - if the file is different to the expected existing format then the conversion will fail and another converter is
        attempted
     - if successful, then the converter function name that was successful is returned
     - if unsuccessful, then 'NONE' is returned

    Liquepy CPT format
    ------------------

    The liquepy format is a simple standard format.

     - File is saved as a CSV
     - all units in metres and MPa.
     - CPT measurements start at line 25
     - Comma separation using ','
     - First 23 lines contain meta data
     - Pre-drill depth is defined in meta data as 'Pre-Drill:,<pre-drill depth>,'
     - Ground water level is defined in meta data as 'Assumed GWL:,<ground water depth>,'
     - Cone area ratio is defined in meta data as 'aratio,<cone area ratio>,'
     - Second line of metadata contains name of converter function (useful for debugging)

    Parameters
    ----------
    :param ffp: str
        Full file path to original CPT file
    :param out_fp: str
        Output folder path where formatted CPT file should be saved
    :param verbose: bool or int
        if true then print to algorithm steps to console
    :return:
    """

    fns = [convert_ags, convert_raw01_xls, convert_raw01_xls_v2, convert_raw01_xls_v3, convert_raw01_xlsx,
           convert_raw02_xlsx_or_raw_03_xls, convert_raw01_xlsx_v2, convert_fugro_raw01, convert_raw01_fugro_v2,
           convert_raw01_extra_cols, convert_tabulated_kpas, convert_raw01_xlsx_from_csv, convert_raw_txt_xlsx,
           convert_raw01_w_kpa, convert_raw01_fugro_xls_v3, convert_raw01_xlsx_space_sep, convert_shortened_ags,
           convert_raw01_w_underscores, convert_raw01_w_meta_on_right]
    for i, fn in enumerate(fns):
        if verbose:
            print(fn.__name__)
        res = fn(ffp, out_fp, verbose=verbose)
        if res == 1:
            return fn.__name__
        if i == len(fns) - 1:
            if verbose:
                print("NOT CONVERTED!")
            return "NONE"


def convert_folder(in_fp, out_fp, verbose=0):
    """
    Reads through a folder of NZGD CPT files and converts them to liquepy format

    :param in_fp: str
        Folder path of existing NZGD CPT files
    :param out_fp: str
        Folder path where formatted NZGD files should be saved
    :param verbose: int or bool
        if true then print to algorithm steps to console
    :return:
    """
    results = []
    ffps = glob.glob(in_fp + "*.xls*")  # Does not handle .csv and .txt
    ffps.sort()

    if not os.path.exists(out_fp):
        os.makedirs(out_fp)
    for ffp in ffps:
        if verbose:
            print(ffp)

        converter_name = convert_file(ffp, out_fp, verbose)
        fname = ffp.split(in_fp)[-1]

        results.append("{0},{1}".format(fname, converter_name))
    if verbose:
        print('SUMMARY: ')
        for line in results:
            print(line)
    ofile = open(out_fp + 'last_processed.txt', 'w')
    ofile.write("\n".join(results))
    ofile.close()
    return results


def trim_missing_at_end_data_df(df_data, neg_lim=None):
    """
    Removes rows at end of file that have empty data

    :param df_data:
    :param neg_lim:
    :return:
    """
    df_data = df_data.reset_index(drop=True)
    nan_rows = df_data[df_data.iloc[:, :3].isnull().T.any().T]
    nan_indexes = list(nan_rows.index)
    if len(nan_indexes):
        i_start = None
        i_end = None
        if nan_indexes[0] == 0:
            i_start = nan_indexes[-1] + 1
            for i in range(1, len(nan_indexes)):
                if nan_indexes[i] - nan_indexes[i - 1] > 1:
                    i_start = nan_indexes[i - 1] + 1
                    i_end = nan_indexes[i]
                    break
        else:
            i_end = nan_indexes[0]
        df_data = df_data[i_start:i_end]

    if neg_lim is not None:
        # remove large neg values
        neg_rows = df_data[df_data.iloc[:, 2] < neg_lim]
        neg_indexes = list(neg_rows.index)
        if len(neg_indexes):
            df_data = df_data[:neg_indexes[0]]

    return df_data


def clean_top(df_top):
    """
    Removes empty data at top of file

    :param df_top:
    :return:
    """
    for i in range(len(df_top)):
        for j in range(4):
            if isinstance(df_top.iloc[i, j], str) and '\r\n' in df_top.iloc[i, j]:
                df_top.iloc[i, j] = df_top.iloc[i, j].replace('\r\n', '')


def convert_ags(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_AGS" in ffp:
        cpt_num = int(cpt_num.split("_AGS")[0])
    elif "_Raw" in ffp:  # e.g. CPT_27240_Raw01.xls
        cpt_num = int(cpt_num.split("_Raw")[0])
    elif "McM" in ffp or "M" in ffp or "DAL" in ffp or "Fitz" in ffp or "ANZ" in ffp:
        cpt_num = cpt_num[:-4]
    else:
        try:
            if cpt_num[-5:] == '.xlsx':
                cpt_num = int(cpt_num[:-5])
            elif cpt_num[-4:] == '.xls':
                cpt_num = int(cpt_num[:-4])
            else:
                return 0
        except ValueError:
            return 0

    df = pd.read_excel(ffp, sheet_name=0)
    cols = list(df)
    if len(df) < 20 or len(cols) < 4:
        return 0
    if df.iloc[20, 0] == D_STR and df.iloc[20, 1] == QC_STR and df.iloc[20, 2] == FS_STR and df.iloc[20, 3] == U2_STR:
        pass
    else:
        return 0
    df_titles = df.iloc[20:21, 0:4]
    df_data = df.iloc[21:, 0:4]
    df_data = trim_missing_at_end_data_df(df_data, neg_lim=-100)

    df_top = df.iloc[0:19, 0:4]
    more = 3
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    clean_top(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]

    df_new = pd.concat([df_top, df_z, df_titles, df_data])
    # df_new.reindex(index=np.arange(len(df_new)))
    # df_new = pd.concat([df_top, df_z, df_data])
    df_top.iloc[3, 1] = 'H'
    df_tt = pd.concat([df_top, df_z, df_titles])

    df_new.to_csv(out_fp + "CPT_{0}.csv".format(cpt_num), index=False)
    return 1


def convert_shortened_ags(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:  # e.g. CPT_27240_Raw01.xls
        cpt_num = int(cpt_num.split("_Raw")[0])
    elif "McM" in ffp or "M" in ffp or "DAL" in ffp or "Fitz" in ffp or "ANZ" in ffp:
        cpt_num = cpt_num[:-4]
    else:
        try:
            if cpt_num[-5:] == '.xlsx':
                cpt_num = int(cpt_num[:-5])
            elif cpt_num[-4:] == '.xls':
                cpt_num = int(cpt_num[:-4])
            else:
                return 0
        except ValueError:
            return 0

    df = pd.read_excel(ffp, sheet_name=0)
    cols = list(df)
    if len(df) < 20 or len(cols) < 4:
        return 0
    ln0 = 14
    ln1 = 15
    if df.iloc[ln0, 0] == 'Depth [m]' and df.iloc[ln0, 1] == 'qc [MPa]' and df.iloc[ln0, 2] == 'fs [MPa]' and df.iloc[ln0, 3] == 'u2 [MPa]':
        ln = ln0
    elif df.iloc[ln1, 0] == 'Depth [m]' and df.iloc[ln1, 1] == 'qc [MPa]' and df.iloc[ln1, 2] == 'fs [MPa]' and df.iloc[ln1, 3] == 'u2 [MPa]':
        ln = ln1
    else:
        return 0

    df_data = df.iloc[ln + 1:, 0:4]
    df_data = trim_missing_at_end_data_df(df_data, neg_lim=-100)

    df_top = df.iloc[0:ln - 1, 0:4]
    heads = [D_STR, QC_STR, FS_STR, U2_STR]
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    limit = 22
    more = limit - len(df_top)
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Predrill":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])
    # df_new = pd.concat([df_top, df_z, df_data])

    df_new.to_csv(out_fp + "CPT_{0}.csv".format(cpt_num), index=False)
    return 1


def convert_tabulated_kpas(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:  # e.g. CPT_27240_Raw01.xls
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    df = pd.read_excel(ffp, sheet_name=0)
    cols = list(df)
    if len(df) < 21 or len(cols) < 4:
        return 0
    if df.iloc[21, 0] == D_STR and df.iloc[21, 1] == QC_STR and df.iloc[21, 2] == 'fs (kPa)' and df.iloc[21, 3] == 'u (kPa)':
        pass
    else:
        return 0
    # df_titles = df.iloc[20:21, 0:4]
    df_data = df.iloc[22:, 0:4]
    df_data = trim_missing_at_end_data_df(df_data, neg_lim=-100)
    df_data.iloc[:, 2] = df_data.iloc[:, 2] / 1e3  # convert to MPa
    df_data.iloc[:, 3] = df_data.iloc[:, 3] / 1e3  # convert to MPa

    df_top = df.iloc[0:20, 0:4]

    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    more = 2
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    for i in range(len(df_top)):
        for j in range(3):
            if df_top.iloc[i, j] == 'Nil':
                df_top.iloc[i, j] = 0.0
            if df_top.iloc[i, j] == 'nil':
                df_top.iloc[i, j] = 0.0
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])
    # df_new = pd.concat([df_top, df_z, df_data])

    df_new.to_csv(out_fp + "CPT_{0}.csv".format(cpt_num), index=False)
    return 1


def convert_raw01_extra_cols(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])

    d_str = 'Penetration Depth (m)'
    # try info all on one page
    try:
        df = pd.read_excel(ffp, sheet_name=0)
        cols = list(df)
        if len(df) < 20 or len(cols) < 4:
            raise ValueError
        if df.iloc[20, 1] == d_str and df.iloc[20, 2] == QC_STR and df.iloc[20, 3] == FS_STR and df.iloc[20, 4] == U2_STR:
            df_data = df.iloc[21:, 1:5]
            df_top = df.iloc[0:19, 0:4]
        else:
            raise ValueError
    except ValueError:
        # try info split over two sheets
        try:
            df = pd.read_excel(ffp, sheet_name=1)
        except IndexError:
            return 0
        cols = list(df)
        if cols[0] == d_str and cols[1] == QC_STR and cols[2] == FS_STR and cols[3] == U2_STR:
            df_data = df.iloc[:, :4]
            df_top = pd.read_excel(ffp, sheet_name=0)
            df_top = df_top.iloc[:min(22, len(df_top)), :4]
        else:
            return 0

    heads = [D_STR, QC_STR, FS_STR, U2_STR]

    df_data = trim_missing_at_end_data_df(df_data)

    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    more = 22 - len(df_top)
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data], sort=True)

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_fugro_raw01(ffp, out_fp, verbose=0):  # e.g. CPT_27643_Raw01.xlsx
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    df = pd.read_excel(ffp, sheet_name=0)
    cols = list(df)
    tline = 44
    uline = 45
    if len(df) < 20 or len(cols) < 4:
        return 0
    if df.iloc[tline, 0] == 'Depth' and df.iloc[tline, 1] == 'Cone' and df.iloc[tline, 2] == 'Friction' and df.iloc[tline, 3] == 'Pore 2':
        if df.iloc[uline, 0] == 'm' and df.iloc[uline, 1] == 'MPa' and df.iloc[uline, 2] == 'MPa' and df.iloc[uline, 3] == 'MPa':
            pass
        else:
            return 0
    else:
        return 0
    df_data = df.iloc[uline + 1:, 0:4]
    df_data = trim_missing_at_end_data_df(df_data)
    df_top = df.iloc[0:19, 0:4]
    more = 3
    zeros = [["", "", "", ""]] * more
    zeros.append([D_STR, QC_STR, FS_STR, U2_STR])
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_data], sort=True)

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xls(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)
    if 'Header' in xf.sheet_names and 'CPT1' in xf.sheet_names:
        if verbose:
            print('found header match at convert_raw01_xls')
    else:
        return 0
    df = pd.read_excel(ffp, sheet_name='CPT1')
    if df.iloc[5, 0] == D_STR and df.iloc[5, 1] == QC_STR and df.iloc[5, 2] == 'fs(kPa)' and df.iloc[5, 3] == 'u (kPa)':
        pass
    else:
        return 0
    df_data = df.iloc[6:, 0:4]
    df_data = trim_missing_at_end_data_df(df_data)
    # Convert columns from kPa to MPa
    df_data.iloc[:, 2] = df_data.iloc[:, 2] / 1e3
    df_data.iloc[:, 3] = df_data.iloc[:, 3] / 1e3

    df_top = pd.read_excel(ffp, sheet_name='Header')
    df_top = df_top.iloc[:, 3:7]

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, 4):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    limit = 20
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top[limit - 2] = df_top[gwl_row]
    if aratio_row is not None and aratio_row > limit:
        df_top[limit - 1] = df_top[aratio_row]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xls_v2(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'Depth (m)'
        qc_str = 'qc Clean (MPa)'
        fs_str = 'fs Clean (MPa)'
        u2_str = 'u2 Clean (MPa)'
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str and cols[3] == u2_str:
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, 0:4]
    # Columns already in MPa
    df_data = trim_missing_at_end_data_df(df_data)

    df_top = pd.read_excel(ffp, sheet_name='CPT_Details')

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, 4):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water Level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    limit = 22
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    # for i in range(more):
    #     df_top.loc[-1] = [0, 0, 0, 0]
    #     df_top.append([0, 0, 0, 0])
    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_w_meta_on_right(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'Depth (m)'
        qc_str = 'qcClean (MPa)'
        fs_str = 'fsClean (MPa)'
        u2_str = 'u2Clean (MPa)'
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str and cols[3] == u2_str:
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, 0:4]
    # Columns already in MPa
    df_data = trim_missing_at_end_data_df(df_data)

    df_top = df.iloc[:, 5:]
    r1 = []
    r2 = []
    r3 = []
    for head in df_top.columns:
        r1.append([head, df_top[head].iloc[0]])
        if not pd.isnull(df_top[head].iloc[2]):
            r2.append([df_top[head].iloc[2], df_top[head].iloc[3]])
        if not pd.isnull(df_top[head].iloc[6]):
            r3.append([df_top[head].iloc[6], df_top[head].iloc[7]])
    dd = r1 + r2 + r3
    df_top = pd.DataFrame(dd)

    n_header_cols = 2
    if n_header_cols < 4:
        for i in range(n_header_cols, 4):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water Level":  # not used
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "ConeAreaRatio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    limit = 22
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))

    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_fugro_xls_v3(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'Depth'
        qc_str = 'Cone'
        fs_str = 'Friction'
        u2_strs = ['Pore 2', 'Pore2']
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str and cols[3] in u2_strs:
            if df.iloc[0, 0] == 'm' and df.iloc[0, 1] == 'MPa' and df.iloc[0, 2] == 'MPa' and df.iloc[0, 3] == 'MPa':
                found = name
            else:
                return 0
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[1:, 0:4]
    # Columns already in MPa
    df_data = trim_missing_at_end_data_df(df_data)

    df_top = pd.read_excel(ffp, sheet_name='Description')

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, 4):
            df_top["C%i" % i] = ""
    else:
        df_top = df_top.iloc[:, :4]
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Derived GWL (m):":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Tip Net Area Ratio:":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
        if df_top.iloc[i, 0] == "Pre-Drill (m):":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            predrill_row = i
    limit = 22
    more = limit - pp

    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_fugro_v2(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    opt0 = ['Depth [m]', 'qc [MPa]', 'fs [MPa]', 'u2 [MPa]']
    opt1 = ['Depth [m]', 'Qc [MPa]', 'Fs [MPa]', 'U2 [MPa]']
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == opt0[0] and cols[1] == opt0[1] and cols[2] == opt0[2] and cols[3] == opt0[3]:
            found = name
        if cols[0] == opt1[0] and cols[1] == opt1[1] and cols[2] == opt1[2] and cols[3] == opt1[3]:
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, 0:4]
    # Convert columns from kPa to MPa
    df_data.iloc[:, 2] = df_data.iloc[:, 2]
    df_data.iloc[:, 3] = df_data.iloc[:, 3]
    df_data = trim_missing_at_end_data_df(df_data)

    if 'TESTING SUMMARY SHEET' in sheet_names:
        df_top = pd.read_excel(ffp, sheet_name='TESTING SUMMARY SHEET')
    elif 'Description' in sheet_names:
        df_top = pd.read_excel(ffp, sheet_name='Description')
    elif 'CPT SUMMARY SHEET' in sheet_names:
        df_top = pd.read_excel(ffp, sheet_name='CPT SUMMARY SHEET')
        dd = []
        for head in df_top.columns:
            dd.append([head, df_top[head].iloc[0]])
        df_top = pd.DataFrame(dd)
    else:
        raise ValueError("Sheet name not detected")

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, 4):
            df_top["C%i" % i] = ""
    elif n_cols > 4:
        df_top = df_top.iloc[:, :4]
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water Level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    limit = 22
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [["", "", "", ""]] * more

    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_w_kpa(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    opt1 = ['Depth [m]', 'Qc [MPa]', 'Fs [KPa]', 'U2 [KPa]']
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        if len(df) < 32:
            continue
        cols = list(df)

        if len(cols) < 3:
            continue
        if df.iloc[31, 0] == opt1[0] and df.iloc[31, 1] == opt1[1] and df.iloc[31, 2] == opt1[2] and df.iloc[31, 3] == opt1[3]:
            found = name
            break
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df_data = df.iloc[32:, 0:4]
    df_data.iloc[:, 2] = df_data.iloc[:, 2] / 1e3
    df_data.iloc[:, 3] = df_data.iloc[:, 3] / 1e3
    df_top = df.iloc[0:19, 0:4]

    heads = [D_STR, QC_STR, FS_STR, U2_STR]

    df_data = trim_missing_at_end_data_df(df_data)

    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    more = 22 - len(df_top)
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_top.index = df_top.index + 1  # shifting index
    df_top = df_top.sort_index()  # sorting by index
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data], sort=True)

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xls_v3(ffp, out_fp, verbose=0):
    d_str = 'Test Length (m)'
    qc_str = 'Cone Resistance (MPa)'
    fs_str = 'Local Friction (MPa)'
    u2_str = 'Pore Pressure (MPa)'

    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 3
    offset = None
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)

        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str:
            if cols[3] == u2_str:
                n_data_cols = 4
            offset = 0
            found = name
        elif cols[1] == d_str and cols[2] == qc_str and cols[3] == fs_str:
            if cols[4] == u2_str:
                n_data_cols = 4
            offset = 1
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, offset:n_data_cols + offset]
    df_data = trim_missing_at_end_data_df(df_data)
    df_data.fillna(0, inplace=True)
    # Convert columns from kPa to MPa
    df_data.iloc[:, 2] = df_data.iloc[:, 2]
    if n_data_cols == 4:
        df_data.iloc[:, 3] = df_data.iloc[:, 3]

    df_top = pd.read_excel(ffp, sheet_name='General')
    dd = []
    for head in df_top.columns:
        dd.append([head, df_top[head].iloc[0]])
    df_top = pd.DataFrame(dd)

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, n_data_cols):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water Level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
        if df_top.iloc[i, 0] == "Predrilled":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            predrill_row = i
    limit = 21
    more = limit - pp

    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    heads = [D_STR, QC_STR, FS_STR, U2_STR]
    if n_data_cols == 3:
        heads = [D_STR, QC_STR, FS_STR]
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_w_underscores(ffp, out_fp, verbose=0):
    d_str = 'Depth_m'
    qc_str = 'qc_mpa'
    fs_str = 'fs_mpa'
    u2_str = 'u2_mpa'

    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 3
    offset = None
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)

        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str:
            if cols[3] == u2_str:
                n_data_cols = 4
            offset = 0
            found = name
        elif cols[1] == d_str and cols[2] == qc_str and cols[3] == fs_str:
            if cols[4] == u2_str:
                n_data_cols = 4
            offset = 1
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, offset:n_data_cols + offset]
    df_data = trim_missing_at_end_data_df(df_data)
    df_data.fillna(0, inplace=True)

    df_top = pd.read_excel(ffp, sheet_name='CPT_Setup')
    dd = []
    for head in df_top.columns:
        dd.append([head, df_top[head].iloc[0]])
    df_top = pd.DataFrame(dd)

    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, n_data_cols):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "WaterLevel":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "ConeAreaRatio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
        if df_top.iloc[i, 0] == "Predrill":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            if df_top.iloc[i, 1] == '' or pd.isnull(df_top.iloc[i, 1]):
                df_top.iloc[i, 1] = '0.0'
            predrill_row = i
    limit = 21
    more = limit - pp

    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    heads = [D_STR, QC_STR, FS_STR, U2_STR]
    if n_data_cols == 3:
        heads = [D_STR, QC_STR, FS_STR]
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xlsx_v2(ffp, out_fp, verbose=0):  # e.g. CPT_38858_Raw01.xlsx
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 3
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'Depth'
        qc_str = 'qc'
        fs_str = 'fs'
        u2_str = 'u2'
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str:
            if df.iloc[0, 1] != '[MPa]' and df.iloc[0, 2] != '[MPa]':
                return 0
            if cols[3] == u2_str:
                n_data_cols = 4
                if df.iloc[0, 3] != '[MPa]':
                    return 0
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[2:, 0:n_data_cols]
    df_data.fillna(0, inplace=True)
    # Convert columns from kPa to MPa
    df_data.iloc[:, 2] = df_data.iloc[:, 2]
    if n_data_cols == 4:
        df_data.iloc[:, 3] = df_data.iloc[:, 3]

    df_data = trim_missing_at_end_data_df(df_data)
    df_top = pd.read_excel(ffp, sheet_name='Header')
    dd = []
    for head in df_top.columns:
        dd.append([head, df_top[head].iloc[0]])
    df_top = pd.DataFrame(dd)

    n_cols = len(list(df_top))
    if n_cols < n_data_cols:
        for i in range(n_cols, n_data_cols):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Waterlevel:":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    limit = 22
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    if n_data_cols == 4:
        heads = [D_STR, QC_STR, FS_STR, U2_STR]
    else:
        heads = [D_STR, QC_STR, FS_STR]
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xlsx_space_sep(ffp, out_fp, verbose=0):  # e.g. CPT_38858_Raw01.xlsx
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 3
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        cols = list(df)
        if len(cols) != 1:
            continue
        titles = df.iloc[9, 0]
        units = df.iloc[10, 0]
        if 'Depth     qc' in titles:
            if '[MPa]     [MPa]' in units:
                found = name
                break

    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    data = []
    for i in range(11, len(df)):
        data.append(df.iloc[i, 0].split()[:4])

    top = []
    for i in range(9):
        top.append(df.iloc[i, 0].split(':'))

    more = 13
    zeros = [[""] * n_data_cols] * more
    heads = [D_STR, QC_STR, FS_STR, U2_STR]
    head_lines = [heads]
    dd = top + zeros + head_lines + data

    df_new = pd.DataFrame(dd)
    df_new.iloc[0, -1] = inspect.stack()[0][3]

    for i in range(len(df_new)):
        if df_new.iloc[i, 0] == "Waterlevel":
            df_new.iloc[i, 0] = 'Assumed GWL:'
        if df_new.iloc[i, 0] == "Cone Area Ratio":
            df_new.iloc[i, 0] = 'aratio'
        if df_new.iloc[i, 0] == "PreDrill":
            df_new.iloc[i, 0] = 'Pre-Drill:'

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw_txt_xlsx(ffp, out_fp, verbose=0):  # e.g. CPT_49062
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 4
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'Depth'
        qc_str = 'qc'
        fs_str = 'fs'
        u2_str = 'u2'
        cols = list(df)
        if len(cols) < 3:
            continue
        if len(df) <= 12:
            continue

        if df.iloc[12, 0] == d_str and df.iloc[12, 1] == qc_str and df.iloc[12, 2] == fs_str and df.iloc[12, 3] == u2_str:
            found = name
        else:
            return 0

    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[15:, 0:n_data_cols]

    df_data = trim_missing_at_end_data_df(df_data)
    df_data.fillna(0, inplace=True)
    df_top = df.iloc[:12, 0:n_data_cols]

    n_cols = len(list(df_top))
    if n_cols < n_data_cols:
        for i in range(n_cols, n_data_cols):
            df_top["C%i" % i] = ""
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Waterlevel:":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
        if df_top.iloc[i, 0] == "PreDrill":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            predrill_row = i
    limit = 22
    more = limit - pp

    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    if n_data_cols == 4:
        heads = [D_STR, QC_STR, FS_STR, U2_STR]
    else:
        heads = [D_STR, QC_STR, FS_STR]
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xlsx_from_csv(ffp, out_fp, verbose=0):  # e.g. CPT_29585_Raw01.xlsx
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 4
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)
        d_str = 'H [m]'
        qc_str = 'qc [MPa]'
        fs_str = 'fs [MPa]'
        u2_str = 'u2 [MPa]'
        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str:
            found = name
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    df_data = df.iloc[:, 0:n_data_cols]
    df_data.fillna(0, inplace=True)

    df_data = trim_missing_at_end_data_df(df_data)

    more = 22

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list([""] * n_data_cols))
    if n_data_cols == 4:
        heads = [D_STR, QC_STR, FS_STR, U2_STR]
    else:
        heads = [D_STR, QC_STR, FS_STR]
    df_headers = pd.DataFrame([heads], columns=[""] * n_data_cols)
    df_data.columns = [""] * n_data_cols
    df_z.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw02_xlsx_or_raw_03_xls(ffp, out_fp, verbose=0):
    cpt_num = ffp.split("CPT_")[-1]
    if "_Raw" in ffp:
        cpt_num = int(cpt_num.split("_Raw")[0])
    else:
        return 0

    xf = pd.ExcelFile(ffp)

    sheet_names = xf.sheet_names
    found = 0
    n_data_cols = 3
    opt1 = ['Test length[1]', 'Cone resistance[2]', 'Local friction[3]', 'Pore pressure u2[6]']
    opt2 = ['Test length[1] m', 'Cone resistance[2] MPa', 'Local friction[3] MPa', 'Pore pressure u2[6] MPa']
    units_in_title = 0
    for name in sheet_names:
        df = pd.read_excel(ffp, sheet_name=name)

        cols = list(df)
        if len(cols) < 3:
            continue
        if cols[0] == opt1[0] and cols[1] == opt1[1] and cols[2] == opt1[2]:
            if df.iloc[0, 0] == 'm' and df.iloc[0, 1] == 'MPa' and df.iloc[0, 2] == 'MPa':
                if cols[3] == opt1[3]:
                    n_data_cols = 4
                found = name
                break
        elif cols[0] == opt2[0] and cols[1] == opt2[1] and cols[2] == opt2[2] and cols[3] == opt2[3]:
            found = name
            units_in_title = 1
            break
    if not found:
        if verbose:
            print("CPT headers not found")
        return 0
    df = pd.read_excel(ffp, sheet_name=found)
    if units_in_title:
        df_data = df.iloc[:, 0:n_data_cols]
    else:
        df_data = df.iloc[1:, 0:n_data_cols]
    df_data.fillna(0, inplace=True)

    df_data = trim_missing_at_end_data_df(df_data)
    if 'General' in sheet_names:
        df_top = pd.read_excel(ffp, sheet_name='General')
        dd = []
        for head in df_top.columns:
            dd.append([head, df_top[head].iloc[0]])
        df_top = pd.DataFrame(dd)
    elif 'Header' in sheet_names:
        df_top = pd.read_excel(ffp, sheet_name='Header')
    else:
        return 0

    n_cols = len(list(df_top))
    if n_cols < n_data_cols:
        for i in range(n_cols, n_data_cols):
            df_top["C%i" % i] = ""
    if n_cols > n_data_cols:
        df_top = df_top.iloc[:, :n_data_cols]
    pp = len(df_top)

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Water Level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Cone Area Ratio":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
        if df_top.iloc[i, 0] == "Predrilled":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            predrill_row = i
    limit = 22
    more = limit - pp

    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]

    zeros = [[""] * n_data_cols] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))

    if n_data_cols == 4:
        heads = [D_STR, QC_STR, FS_STR, U2_STR]
        df_data.columns = list(df_top)
    else:
        heads = [D_STR, QC_STR, FS_STR]
        df_data.columns = list(df_top)
    df_headers = pd.DataFrame([heads], columns=list(df_top))
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)
    return 1


def convert_raw01_xlsx(ffp, out_fp, verbose=0):
    if "_Raw" not in ffp:
        return 0
    cpt_num = ffp.split("CPT_")[-1]
    cpt_num = int(cpt_num.split("_Raw")[0])
    xf = pd.ExcelFile(ffp)
    if 'Header' in xf.sheet_names and 'Data' in xf.sheet_names:
        if verbose:
            print('found sheet names at convert_raw01_xlsx')
    else:
        return 0
    df = pd.read_excel(ffp, sheet_name='Data')
    d_str = 'Depth [m]'
    qc_str = 'Cone resistance (qc) in MPa'
    fs_str = 'Sleeve friction (fs) in MPa'
    u2_str = 'Dynamic pore pressure (u2) in MPa'
    cols = list(df)
    if cols[0] == d_str and cols[1] == qc_str and cols[2] == fs_str and cols[3] == u2_str:
        pass
    else:
        return 0
    df_data = df.iloc[:, 0:4]

    df_data = trim_missing_at_end_data_df(df_data)

    df_top = pd.read_excel(ffp, sheet_name='Header')
    df_top = df_top.iloc[:, 0:4]

    pp = len(df_top)
    limit = 22
    more = limit - pp

    gwl_row = None
    aratio_row = None
    predrill_row = None
    for i in range(len(df_top)):
        if df_top.iloc[i, 0] == "Predrill":
            df_top.iloc[i, 0] = 'Pre-Drill:'
            predrill_row = i
        if df_top.iloc[i, 0] == "Water level":
            df_top.iloc[i, 0] = 'Assumed GWL:'
            if df_top.iloc[i, 1] < 0:
                df_top.iloc[i, 1] *= -1
            gwl_row = i
        if df_top.iloc[i, 0] == "Alpha Factor":
            df_top.iloc[i, 0] = 'aratio'
            aratio_row = i
    if predrill_row is not None and predrill_row > limit:
        df_top.iloc[limit - 3, :] = df_top.iloc[predrill_row, :]
    if gwl_row is not None and gwl_row > limit:
        df_top.iloc[limit - 2, :] = df_top.iloc[gwl_row, :]
    if aratio_row is not None and aratio_row > limit:
        df_top.iloc[limit - 1, :] = df_top.iloc[aratio_row, :]
    if more < 0:
        df_top = df_top[:limit]
    n_cols = len(list(df_top))
    if n_cols < 4:
        for i in range(n_cols, 4):
            df_top["C%i" % i] = ""

    n_cols = len(list(df_top))
    if more < 0:
        df_top = df_top[:limit]
    zeros = [["", "", "", ""]] * more
    df_z = pd.DataFrame(zeros, columns=list(df_top))
    df_headers = pd.DataFrame([[D_STR, QC_STR, FS_STR, U2_STR]], columns=list(df_top))
    df_data.columns = list(df_top)
    df_top.iloc[0, -1] = inspect.stack()[0][3]
    df_new = pd.concat([df_top, df_z, df_headers, df_data])

    df_new.to_csv(out_fp + "CPT_%i.csv" % cpt_num, index=False)

    return 1

