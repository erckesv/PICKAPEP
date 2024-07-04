"""PICKAPEP peptide generator

This application allows the user to generate peptides with various modifications and cyclizations.
It also has an option for multiple peptide processing.

"""

import sys
import os
import openpyxl
import pandas as pd
import numpy as np
from os import path
import re
import pickle
import random
import tkinter as tk
import pyperclip
import subprocess
from tkinter import ttk, filedialog, messagebox
from PIL import ImageTk, Image, ImageDraw, ImageFont
from io import BytesIO
from rdkit import Chem, RDLogger
from rdkit.Chem import (
    rdChemReactions,
    Draw,
    MolStandardize,
    AllChem,
    Descriptors,
    rdBase,
    rdCoordGen,
)
from rdkit.Chem.Draw import IPythonConsole

# Close Pyinstaller Splash
try:
    import pyi_splash

    pyi_splash.update_text("UI Loaded ...")
    pyi_splash.close()
except:
    pass

# Suppress RD Kit Errors
RDLogger.DisableLog("rdApp.*")

## DATAPATHDEFINITION

# os.chdir(sys._MEIPASS)
bundle_dir = os.path.abspath(os.path.dirname(__file__))


path_to_AminoAcids = os.path.join(bundle_dir, "Python_AminoAcids.csv")
path_to_TerminalMod = os.path.join(bundle_dir, "Python_TerminalMod.csv")
path_to_CyclReactions = os.path.join(bundle_dir, "Python_CyclReactions.csv")
path_to_ModReactions = os.path.join(bundle_dir, "Python_ModReactions.csv")
path_to_dict_default = os.path.join(bundle_dir, "dict_default.pkl")
path_to_helpbutton = os.path.join(bundle_dir, "Helpbutton_small.png")
path_to_teriyaki = os.path.join(bundle_dir, "Teriyaki.png")
path_to_teriyaki_easteregg = os.path.join(bundle_dir, "Teriyaki_Easteregg.png")
path_to_infobutton = os.path.join(bundle_dir, "Infobutton_small.png")
path_to_HS_Userinput1_aa = os.path.join(
    bundle_dir, "HelpSheet_Userinput1_aminoacids.png"
)
path_to_HS_Userinput1_mpc = os.path.join(
    bundle_dir, "HelpSheet_Userinput1_multiplepepcalc.png"
)
path_to_HS_Userinput1 = os.path.join(bundle_dir, "HelpSheet_Userinput1.png")
path_to_HS_Userinput2 = os.path.join(bundle_dir, "HelpSheet_Userinput2.png")
path_to_HS_Userinput3 = os.path.join(bundle_dir, "HelpSheet_Userinput3.png")
path_to_HS_Userinput2_cycl = os.path.join(bundle_dir, "HelpSheet_Userinput2_cycl.png")
path_to_HS_Userinput2_mod = os.path.join(bundle_dir, "HelpSheet_Userinput2_mod.png")


class Peptide:
    """
    A class to represent a peptide.

    Attributes
    ----------
    peptide : rdkit.Chem.rdchem.Mol
        RDkit mol object of the peptide with assigned mapnums
    mapnumstore : list
        list to store highest mapnums assigned during peptide assembly
    highestmapnum : int
        current highest mapnum
    df_sequence : pd.DataFrame
        DataFrame with sequence information
    df_cycl : pd.DataFrame
        DataFrame with cyclization information
    df_aamod : pd.DataFrame
        DataFrame with modification information
    cycl_num : int
        number of cyclizations conducted
    aamod_num: int
        number of modifications conducted
    pepinfo: pd.DataFrame
        DataFrame containing calculated peptide parameters
    """

    def __init__(self, peptide, mapnumstore, df_sequence, df_cycl, df_aamod, pepinfo):
        self.peptide = peptide  # RDkit mol object with mapnums assigned
        self.mapnumstore = mapnumstore
        self.highestmapnum = mapnumstore[-1]
        self.df_sequence = df_sequence
        self.df_cycl = df_cycl
        self.df_aamod = df_aamod
        self.cycl_num = len(df_cycl)
        self.aamod_num = len(df_aamod)
        self.pepinfo = pepinfo


def errorupdate(
    window_status,
    calc_type,
    calc_table,
    fastcalc_pepnum,
    errormessage,
    new_window_status,
):
    """Update user (calc type == slow) or calculation table (calc type == fast) in case PICKAPEP encounters a userinput it cannot process.

    Parameters
    ----------
    window_status : int
        Number of current window the user is encountering the error from.
    calc_type: str
        String, either 'slow' or 'fast' depending on processing mode.
    calc_table: pd.DataFrame
        Data frame containing in case of calc_type = fast the imported table with peptides to calculate.
    fastcalc_pepnum:
        Number of the current peptide processed in case of calc_type == fast to update the status line with the encountered error.
    errormessage: str
        String containing the information about the error to generate a problem specific errormessage.
    new_window_status: int
        Number of window the user should skip to after beining informed/updated about the error.

    Returns
    ------
    window_status : int
        Updated number of current window the user is directed to after encountering an error.
    calc_table: pd.DataFrame
        Updated data frame containing a case specific errormessage in the status column.

    """

    if calc_type == "slow":
        window_status = new_window_status
        messagebox.showinfo("Info", errormessage)
    else:
        calc_table.at[fastcalc_pepnum - 1, "Status"] = errormessage
        calc_table.iloc[fastcalc_pepnum - 1, 7:] = ""
    return window_status, calc_table


def window_size_centered(window, width, height):
    """Placement of tkinter GUI windows in the center of the screen.

    Parameters
    ----------
    window: tkinter.Tk object
        Tkinter window which has to be centred.
    width: int
        Width of the GUI window.
    height: int
        Height of the GUI window.

    """

    window.geometry("{}x{}".format(width, height))
    window.update_idletasks()
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    x = (screen_width - width) // 2
    y = (screen_height - height) // 2 - 40
    if screen_height < (y + height + 40):
        y = 5
        x = 5
    window.geometry("+{}+{}".format(x, y))
    window.lift()


def create_peptide_summary(peptide_memory):
    """Summarize parameters of the current peptide in a data frame.

    Parameters
    ----------
    peptide_memory: list
        List containing all peptide objects from the linear sequence to subsequent modified/cyclized forms.

    Returns
    ------
    df_summary_properties: pd.DataFrame
        Data frame containing all parameters calculated or summarized for the current/final peptide.

    """

    df_sequence_overview = peptide_memory[-1].df_sequence.copy()
    aa_list = df_sequence_overview["AminoAcidInfo"].tolist()
    cycl_info_list, mod_info_list = get_info_lists(peptide_memory)
    df_summary_properties = peptide_memory[-1].pepinfo
    nterm = peptide_memory[-1].df_sequence["AminoAcid"][0]
    cterm = peptide_memory[-1].df_sequence["AminoAcid"][len(aa_list) - 1]
    sequence_info = ",".join(aa_list)
    nterm_str = nterm if nterm != "none" else ""
    cterm_str = cterm if cterm != "none" else ""
    linearsequence_info = (
        nterm_str
        + "-".join(list(peptide_memory[-1].df_sequence["AminoAcid"][1:-1]))
        + cterm_str
    )
    cyclization_info = "; ".join(
        cycl_info for cycl_info in cycl_info_list if cycl_info != " "
    )
    modification_info = "; ".join(
        mod_info for mod_info in mod_info_list if mod_info != " "
    )
    summary_properties_additional = {
        "Info": [
            "Sequence",
            "Linear Sequence",
            "Cyclizations",
            "Modifications",
            "Smiles Code",
        ],
        "Value": [
            sequence_info,
            linearsequence_info,
            cyclization_info,
            modification_info,
            Chem.MolToSmiles(peptide_memory[-1].peptide),
        ],
    }
    df_summary_properties_additional = pd.DataFrame(summary_properties_additional)
    df_summary_properties = pd.concat(
        [df_summary_properties_additional, df_summary_properties], ignore_index=True
    )
    df_summary_properties.loc[len(df_summary_properties)] = [None, None]
    df_summary_properties.loc[len(df_summary_properties)] = [
        "Calculation with PICKAPEP (Erckes V. et al. 2024, J. Pept. Sci.)",
        None,
    ]

    return df_summary_properties


def get_info_lists(peptide_memory):
    """Summarize each cyclization and modification reaction in a list for fast overview.

    Parameters
    ----------
    peptide_memory: list
        List containing all peptide objects from the linear sequence to subsequent modified/cyclized forms.

    Returns
    ------
    cycl_info_list: list
        List containing detailed description of conducted cyclization reactions.
    mod_info_list: list
        List containing detailed description of conducted modification reactions.

    """

    df_cyclpep = peptide_memory[-1].df_cycl
    cycl_info_list = []
    for num, row in df_cyclpep.iterrows():
        num1_info = row["AA1"]
        num2_info = row["AA2"]
        row_info = (
            "("
            + str(num + 1)
            + ")"
            + num1_info
            + "-"
            + num2_info
            + ": "
            + row["CyclName"]
            + ", "
            + row["CyclType"]
        )
        if row["CyclType"] == "Bicycle Formation":
            num3_info = row["AA3"]
            row_info = (
                "("
                + str(num + 1)
                + ")"
                + num1_info
                + "-"
                + num2_info
                + "-"
                + num3_info
                + ": "
                + row["CyclName"]
                + ", "
                + row["CyclType"]
            )
        if row["CyclType"] == "Stitched":
            stitch_info = row["AA3"]
            row_info = (
                "("
                + str(num + 1)
                + ")"
                + num1_info
                + "-"
                + stitch_info
                + "(stitch) -"
                + num2_info
                + ": "
                + row["CyclName"]
                + ", "
                + row["CyclType"]
            )
        cycl_info_list.append(row_info)

    df_modpep = peptide_memory[-1].df_aamod
    mod_info_list = []

    for num, row in df_modpep.iterrows():
        num_info = row["AAModif"]
        if row["Additional"] == 0:
            additional_modinfo = ""
        else:
            additional_modinfo = row["Additional"] + ", "
        row_info = (
            "("
            + str(num + 1)
            + ")"
            + num_info
            + ": "
            + additional_modinfo
            + row["ModifType"]
        )
        mod_info_list.append(row_info)
    return cycl_info_list, mod_info_list


def generate_threedim_pymol(peptide_memory):
    """Embed mol threedimensional, save it as sdf file and open it in PyMOl (in case PICKAPEP cannot find a path to pyMOL, the user can choose the path manually).

    Parameters
    ----------
    peptide_memory: list
        List containing all peptide objects from the linear sequence to subsequent modified/cyclized forms.

    """

    open_pymol = True
    try:
        molH = AllChem.AddHs(peptide_memory[-1].peptide)
        ps = AllChem.ETKDGv2()
        AllChem.EmbedMolecule(molH, ps)
        sdf_path = os.path.join(bundle_dir, "molecule.sdf")
        Chem.SDWriter(sdf_path).write(molH, confId=-1)
        atom_positions = molH.GetConformer().GetPositions()

        x_coords = []
        y_coords = []
        z_coords = []

        for pos in atom_positions:
            x, y, z = pos
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)

        if (
            all(x == 0.0 for x in x_coords)
            or all(y == 0.0 for y in y_coords)
            or all(z == 0.0 for z in z_coords)
        ):
            raise Exception

    except Exception as e:
        try:
            proceed_3D = messagebox.askyesno(
                "Error",
                "PICKAPEP cannot generate the 3D molecule as expected. \n \nDo you want to try again with a less optimized option? \nThis might take a few minutes.",
            )
            if proceed_3D:
                molH = AllChem.AddHs(peptide_memory[-1].peptide)
                ps = AllChem.ETKDGv2()
                ps.useRandomCoords = True
                AllChem.EmbedMolecule(molH, ps)
                sdf_path = os.path.join(bundle_dir, "molecule.sdf")
                Chem.SDWriter(sdf_path).write(molH, confId=-1)
                atom_positions = molH.GetConformer().GetPositions()
            else:
                open_pymol = False
        except Exception as e:
            tk.messagebox.showwarning(
                message="Unfortunately PICKAPEP can not generate \na 3D molecule from your peptide"
            )
            open_pymol = False

    pymol_path_filename = os.path.join(bundle_dir, "pymol_path.txt")
    while open_pymol == True:
        try:
            with open(pymol_path_filename, "r") as file:
                pymol_path = file.read()
            if pymol_path:
                # Open Pymol as subprocess to display 3D peptide molecule
                subprocess.run([pymol_path, sdf_path])
                open_pymol = False
            else:
                raise Exception
        except:
            open_pymol = False
            proceed = messagebox.askyesno(
                "Error",
                "PICKAPEP cannot open PyMOL! Please choose a valid data path for the PyMOL application. \n \nDo you want to proceed?",
            )
            if proceed:
                pymol_path_new = filedialog.askopenfilename(
                    title="Please choose the PyMOL application"
                )
                with open(pymol_path_filename, "w") as file:
                    file.write(pymol_path_new)
                open_pymol = True


## USER INPUT 1 ##
## Sequence input, N- and C- Terminus, Possibility for Cyclization and Modification ##


def userinput1(
    nterm_opt, cterm_opt, dict_default, window_status, df_aminoacids, default_sequence, path_to_teriyaki
):

    def userinput1_click_enter_button(event=None):
        nonlocal sequence, cterm_opt, cterm, nterm_opt, nterm, cycl, window_status, default_options, default_pept
        sequence = linseq_entry.get()
        cterm = cterm_combobox.get()
        nterm = nterm_combobox.get()
        cycl = furtherproc_var.get()

        upper_count = 0
        for letter in sequence:
            if letter.isupper():
                upper_count += 1

        if upper_count > 100:
            result = messagebox.askquestion(
                "",
                "Peptide sequence > 100 amino acids! Calculation times may get very long. Do you want to proceed?",
            )
            if result == "no":
                return

        if sequence != "":
            if sequence in default_options:
                default_pept = linseq_entry.get()
                if cterm != "none" or nterm != "none":
                    tk.messagebox.showwarning(
                        message="Please note, that N- or C-terminus can \nnot be modified in default peptides"
                    )
            if cterm in cterm_opt and nterm in nterm_opt:
                window_status = 2
                if cycl == False:
                    window_status = 3
                window.destroy()
            else:
                tk.messagebox.showwarning(message="Please choose a valid terminus")
        else:
            tk.messagebox.showwarning(message="Please enter a peptide sequence")

    def userinput1_click_aabutton(aaoption):
        nonlocal aa_backbone_pattern, df_aminoacids_display

        for num_row, row in df_aminoacids_display.iterrows():
            aa_type = row["Title"]
            aa_type_df = df_aminoacids[df_aminoacids["Type"] == aa_type].reset_index(
                drop=True
            )
            pepmols = [Chem.MolFromSmiles(x) for x in aa_type_df["Smiles"].tolist()]
            for idx, mol in enumerate(pepmols):
                if aa_type_df.loc[idx, "Name"] == 0:
                    label_name = ""
                else:
                    label_name = aa_type_df.loc[idx, "Name"]
                label_letter = aa_type_df.loc[idx, "Entry_Letter"]
                if aa_type_df.loc[idx, "Entry_Code"] == 0:
                    label_code = ""
                else:
                    label_code = aa_type_df.loc[idx, "Entry_Code"]
                mol.SetProp("label", f"{label_name}\n{label_letter} {label_code}")
                mol.SetProp("label", f"{label_letter} {label_code}")
            df_aminoacids_display.loc[num_row, "Mols"] = pepmols

        if len(pepmols) == 0 and aaoption == "Own Amino Acids":
            messagebox.showinfo("Info", "No Own Amino Acids added")
            return

        index_option = df_aminoacids_display.index[
            df_aminoacids_display["Title"] == aaoption
        ]
        pepmols = df_aminoacids_display.loc[index_option, "Mols"].values[0]

        AllChem.Compute2DCoords(aa_backbone_pattern)
        for m in pepmols:
            _ = AllChem.GenerateDepictionMatching2DStructure(m, aa_backbone_pattern)
        img = Draw.MolsToGridImage(
            pepmols,
            molsPerRow=int(np.ceil(np.sqrt(len(pepmols)))),
            maxMols=len(pepmols),
            legends=[x.GetProp("label") for x in pepmols],
            subImgSize=(140, 140),
            returnPNG=True,
            useSVG=False,
        )

        subwindow = tk.Toplevel(window)
        subwindow.title(aaoption + " Amino Acids")

        aaoption_display = tk.Frame(subwindow, bg="#f0f0f0")
        aaoption_display.grid(row=0, column=0, sticky="news")

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow.iconphoto(False, Teriyaki)

        img_img = Image.open(BytesIO(img.data))
        img_tk = ImageTk.PhotoImage(img_img)
        option_display = tk.Label(aaoption_display, image=img_tk)
        option_display.grid(row=0, column=0, sticky="news", padx=5, pady=5)

        subwindow.mainloop()

    def click_addaa():
        nonlocal df_aminoacids
        aa_entrycode = addaa_entrycode_entry.get()
        aa_entrycode = aa_entrycode.replace(" ", "")
        aa_smiles = addaa_smiles_entry.get()
        pattern_entrycode = r"^[A-Z][a-z0-9]*$"
        pattern_smiles_rx = "[NX3][CX4][CX3](=[O])[O]"
        pattern_smiles = "[NX3H,NX3H2][CX4][CX3](=[O])[OH]"
        pattern_smiles_cterm = "[NX3H,NX3H2][CX4][CX3](=[O])[O]"
        pattern_smiles_nterm = "[NX3][CX4][CX3](=[O])[OH]"

        if re.match(pattern_entrycode, aa_entrycode) is None:
            tk.messagebox.showwarning(
                message="Nonvalid Entry Code: "
                + aa_entrycode
                + "\nEntry code must contain capital letter followed by lower case letters and/or numbers (e.g. Xx1)"
            )
        elif (
            aa_entrycode in df_aminoacids.loc[:, "Entry_Code"].tolist()
            or aa_entrycode in df_aminoacids.loc[:, "Entry_Letter"].tolist()
        ):
            tk.messagebox.showwarning(
                message="This Entry Code "
                + aa_entrycode
                + " is already taken! \nPlease choose another entry code."
            )
        else:
            try:
                Chem.SanitizeMol(
                    Chem.MolFromSmiles(aa_smiles),
                    Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
                )
                if (
                    len(
                        list(
                            (Chem.MolFromSmiles(aa_smiles)).GetSubstructMatches(
                                Chem.MolFromSmarts(pattern_smiles_rx)
                            )
                        )
                    )
                    > 1
                ):
                    final_add = False
                    tk.messagebox.showwarning(
                        message="Nonvalid Smiles Code: \nAmino acid smiles may only contain one backbone pattern RN(R)C(R)COOR"
                    )
                elif (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                    Chem.MolFromSmarts(pattern_smiles)
                ):
                    final_add = tk.messagebox.askyesno(
                        "Question", "Add " + aa_entrycode + " to available amino acids?"
                    )
                elif (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                    Chem.MolFromSmarts(pattern_smiles_nterm)
                ):
                    final_add = tk.messagebox.askyesno(
                        "Question",
                        "Add "
                        + aa_entrycode
                        + " to available amino acids? \n\nATTENTION: R2N- pattern, this amino acid can only be used N-terminal",
                    )
                elif (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                    Chem.MolFromSmarts(pattern_smiles_cterm)
                ):
                    final_add = tk.messagebox.askyesno(
                        "Question",
                        "Add "
                        + aa_entrycode
                        + " to available amino acids? \n\nATTENTION: -COOR pattern, R will only be kept when amino acid is used C-terminal",
                    )
                else:
                    final_add = False
                    tk.messagebox.showwarning(
                        message="Nonvalid Smiles Code: \nAmino acid smiles code must contain backbone pattern NH(R)C(R)COOH"
                    )
                if final_add:
                    addaa_row = pd.DataFrame(
                        [
                            {
                                "Type": "Own Amino Acids",
                                "Entry_Letter": aa_entrycode,
                                "Smiles": aa_smiles,
                            }
                        ]
                    )
                    df_aminoacids = pd.concat(
                        [df_aminoacids, addaa_row], ignore_index=True
                    ).fillna(0)
            except:
                tk.messagebox.showwarning(
                    message="Nonvalid Smiles Code: \nAmino acid smiles code must contain backbone pattern NH(R)C(R)COOH"
                )

    def click_startcalculation_button():
        nonlocal calc_type, calc_table, calc_table_titles, calc_table_file_path, window_status, cycl
        calc_table_file_path = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx;*.xls"), ("CSV files", "*.csv")]
        )

        if calc_table_file_path:
            try:
                if calc_table_file_path.endswith(".csv"):
                    calc_table = pd.read_csv(
                        calc_table_file_path, skip_blank_lines=False
                    ).fillna(" ")
                elif calc_table_file_path.endswith(
                    ".xlsx"
                ) or calc_table_file_path.endswith(".xls"):
                    calc_table = pd.read_excel(calc_table_file_path).fillna(" ")

                if len(calc_table) > 100:
                    result = messagebox.askquestion(
                        "",
                        "> 100 peptides to calculate! Calculation times may get very long. Do you want to proceed?",
                    )
                    if result == "no":
                        return
                df_width = len(calc_table_titles)
                if calc_table.shape[1] > df_width:
                    calc_table = calc_table.iloc[:, :df_width]
                elif calc_table.shape[1] < df_width:
                    num_columns_to_add = df_width - calc_table.shape[1]
                    for i in range(num_columns_to_add):
                        calc_table[i + calc_table.shape[1], "New_Column"] = " "
                calc_table.columns = calc_table_titles
                window_status = 2
                calc_type = "fast"
                cycl = True
                window.destroy()
            except:
                tk.messagebox.showwarning(
                    message="Not possible to Import Calculation Table"
                )
                return

    def userinput1_click_importaabutton():
        nonlocal df_aminoacids
        aaimport_file_path = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx;*.xls"), ("CSV Files", "*.csv")]
        )
        allowed = r"[^A-Za-z0-9]"
        pattern_entrycode = r"^[A-Z][a-z0-9]*$"
        pattern_smiles = "[NX3H,NX3H2][CX4][CX3](=[O])[OH]"
        pattern_smiles_cterm = "[NX3H,NX3H2][CX4][CX3](=[O])[O]"
        pattern_smiles_nterm = "[NX3][CX4][CX3](=[O])[OH]"
        redundantaa = ""

        if aaimport_file_path:
            try:
                if aaimport_file_path.endswith(".csv"):
                    df_importfile = pd.read_csv(
                        aaimport_file_path, usecols=[0, 1], header=0
                    )
                elif aaimport_file_path.endswith(
                    ".xlsx"
                ) or aaimport_file_path.endswith(".xls"):
                    df_importfile = pd.read_excel(aaimport_file_path, usecols=[0, 1])
                addaacount = 0
                aatoadd = len(df_importfile)
            except:
                tk.messagebox.showwarning(message="Not possible to Import File")
                return

        for index, row in df_importfile.iterrows():
            aa_entrycode = row.iloc[0]
            aa_entrycode = re.sub(allowed, "", aa_entrycode)
            aa_smiles = row.iloc[1]
            aa_smiles = aa_smiles.replace(" ", "")
            if re.match(pattern_entrycode, aa_entrycode) is None:
                continue
            elif (
                aa_entrycode in df_aminoacids.loc[:, "Entry_Code"].tolist()
                or aa_entrycode in df_aminoacids.loc[:, "Entry_Letter"].tolist()
            ):
                redundantaa += ", " + aa_entrycode
            else:
                try:
                    Chem.SanitizeMol(
                        Chem.MolFromSmiles(aa_smiles),
                        Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
                    )
                    if (
                        (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                            Chem.MolFromSmarts(pattern_smiles)
                        )
                        or (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                            Chem.MolFromSmarts(pattern_smiles_nterm)
                        )
                        or (Chem.MolFromSmiles(aa_smiles)).HasSubstructMatch(
                            Chem.MolFromSmarts(pattern_smiles_cterm)
                        )
                    ):
                        addaa_row = pd.DataFrame(
                            [
                                {
                                    "Type": "Own Amino Acids",
                                    "Entry_Letter": aa_entrycode,
                                    "Smiles": aa_smiles,
                                }
                            ]
                        )
                        df_aminoacids = pd.concat(
                            [df_aminoacids, addaa_row], ignore_index=True
                        ).fillna(0)
                        addaacount += 1
                except:
                    continue
        aaaddinfo = str(addaacount) + "/" + str(aatoadd) + " amino acids added"
        if redundantaa:
            aaaddinfo = (
                str(addaacount)
                + "/"
                + str(aatoadd)
                + " amino acids added \n\nFollowing entry codes already exist:\n"
                + redundantaa[2:]
            )
        messagebox.showinfo("Info", aaaddinfo)

    def userinput1_click_exportaabutton():

        nonlocal df_aminoacids

        df_ownaminoacids = df_aminoacids[
            df_aminoacids["Type"] == "Own Amino Acids"
        ].reset_index(drop=True)
        if len(df_ownaminoacids) == 0:
            messagebox.showinfo("Info", "No Own Amino Acids added")
            return
        df_ownaminoacids = df_ownaminoacids[["Entry_Letter", "Smiles"]]
        save_file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=(
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ),
        )
        if save_file_path:
            if save_file_path.endswith(".csv"):
                df_ownaminoacids.to_csv(save_file_path, index=False)
            elif save_file_path.endswith(".xlsx"):
                df_ownaminoacids.to_excel(save_file_path, index=False)
            messagebox.showinfo("Save", "File saved successfully")

    def help_button1_click(event=None):

        subwindow_help1 = tk.Toplevel(window)
        subwindow_help1.title("PICKAPEP SUPPORT")

        help1_display = tk.Frame(subwindow_help1, bg="#ffffff")
        help1_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help1.iconphoto(False, Teriyaki)

        img_Help1 = Image.open(path_to_HS_Userinput1)
        img_Help1_width, img_Help1_height = img_Help1.size
        ratio_size = img_Help1_height / img_Help1_width
        target_width = 800
        img_Help1 = img_Help1.resize((target_width, int(target_width * ratio_size)))
        img_Help1_tk = ImageTk.PhotoImage(img_Help1)
        Help1_display = tk.Label(help1_display, image=img_Help1_tk, bg="#ffffff")
        Help1_display.pack(padx=50, pady=50)

        subwindow_help1.mainloop()

    def help_button1_aminoacids_click():

        subwindow_help1_aa = tk.Toplevel(window)
        subwindow_help1_aa.title("PICKAPEP SUPPORT")

        help1_aa_display = tk.Frame(subwindow_help1_aa, bg="#ffffff")
        help1_aa_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help1_aa.iconphoto(False, Teriyaki)

        img_Help1_aa = Image.open(path_to_HS_Userinput1_aa)
        img_Help1_aa_width, img_Help1_aa_height = img_Help1_aa.size
        ratio_size = img_Help1_aa_height / img_Help1_aa_width
        target_width = 800
        img_Help1_aa = img_Help1_aa.resize(
            (target_width, int(target_width * ratio_size))
        )
        img_Help1_aa_tk = ImageTk.PhotoImage(img_Help1_aa)
        Help1_aa_display = tk.Label(
            help1_aa_display, image=img_Help1_aa_tk, bg="#ffffff"
        )
        Help1_aa_display.pack(padx=50, pady=50)

        subwindow_help1_aa.mainloop()

    def help_button1_multiplepepcalc_click(event=None):

        subwindow_help1_mpc = tk.Toplevel(window)
        subwindow_help1_mpc.title("PICKAPEP SUPPORT")

        help1_mpc_display = tk.Frame(subwindow_help1_mpc, bg="#ffffff")
        help1_mpc_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help1_mpc.iconphoto(False, Teriyaki)

        img_Help1_mpc = Image.open(path_to_HS_Userinput1_mpc)
        img_Help1_mpc_width, img_Help1_mpc_height = img_Help1_mpc.size
        ratio_size = img_Help1_mpc_height / img_Help1_mpc_width
        target_width = 800
        img_Help1_mpc = img_Help1_mpc.resize(
            (target_width, int(target_width * ratio_size))
        )
        img_Help1_mpc_tk = ImageTk.PhotoImage(img_Help1_mpc)
        Help1_mpc_display = tk.Label(
            help1_mpc_display, image=img_Help1_mpc_tk, bg="#ffffff"
        )
        Help1_mpc_display.pack(padx=50, pady=50)

        subwindow_help1_mpc.mainloop()

    def click_close_button(event=None):
        nonlocal window_status
        done_input = messagebox.askquestion(
            " ",
            "Are you sure you want to close the program? \n \nATTENTION: Self-added amino acids will be lost! \nSave custom amino acids before closing the program.",
        )
        if done_input == "yes":
            window_status = 4
            window.destroy()

    sequence = cterm = nterm = ""
    cycl = False
    calc_type = "slow"
    calc_table = pd.DataFrame()
    calc_table_titles = [
        "Name",
        "Input Linear Sequence",
        "Input Nterm",
        "Input Cterm",
        "Input Cyclizations",
        "Input Modifications",
        "Status",
        "Sequence",
        "Linear Sequence",
        "Cyclizations",
        "Modifications",
        "Smiles Code",
        "Chemical Formula",
        "Molecular Weight",
        "Exact Weight",
        "[M+H]+",
        "[M+2H]2+",
        "[M+3H]3+",
        "[M+4H]4+",
        "Number of Amino Acids",
        "LogP",
        "TPSA",
        "Rotatable Bond",
        "Number of HBA",
        "Number of HBD",
        "Aromatic Rings",
    ]
    calc_table_file_path = ""
    aa_backbone_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[O])[O]")

    window = tk.Tk()
    window.title("PICKAPEP")
    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    info_png = Image.open(path_to_infobutton)
    resize_info_png = info_png.resize((20, 20))
    info_button_png = ImageTk.PhotoImage(resize_info_png)

    window_frame = tk.Frame(window, bg="#f0f0f0")
    window_frame.pack()

    window_frame.columnconfigure(0, minsize=500)
    window_frame.columnconfigure(1, minsize=500)

    entry_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    entry_frame.grid(row=0, column=0, sticky="news", padx=50, pady=100)

    linearpeptidetitle_label = tk.Label(
        entry_frame, text="LINEAR PEPTIDE ASSEMBLY", bg="#d9d9d9"
    )
    linearpeptidetitle_label.pack(side=tk.TOP, fill=tk.BOTH, padx=30, pady=(30, 0))

    linseq_frame = tk.LabelFrame(entry_frame, text="")
    linseq_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=30, pady=(15, 15))

    linseq_label = tk.Label(linseq_frame, text="Peptide Sequence")
    linseq_label.pack(pady=(10, 0))

    default_pept = ""
    default_options = list(dict_default.keys())
    linseq_entry = ttk.Combobox(linseq_frame, values=default_options, width=50)
    linseq_entry.set(default_sequence)
    linseq_entry.pack(padx=15, pady=(10, 15))

    termmod_frame = tk.LabelFrame(entry_frame, text="")
    termmod_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=30, pady=15)

    termmod_label = tk.Label(termmod_frame, text="Choose N- and C- Terminus")
    termmod_label.grid(row=0, column=0, columnspan=2, pady=10)

    nterm_label = tk.Label(termmod_frame, text="N-Terminus (Tail)")
    nterm_label.grid(row=1, column=0, pady=5)

    nterm_combobox = ttk.Combobox(termmod_frame, values=nterm_opt)
    nterm_combobox.grid(row=1, column=1, pady=5)
    nterm_combobox.set("none")

    cterm_label = tk.Label(termmod_frame, text="C-Terminus (Head)")
    cterm_label.grid(row=2, column=0, pady=(5, 15))

    cterm_combobox = ttk.Combobox(termmod_frame, values=cterm_opt)
    cterm_combobox.grid(row=2, column=1, pady=(5, 15))
    cterm_combobox.set("none")

    for widget in termmod_frame.winfo_children():
        widget.grid_configure(padx=20)

    furtherproc_frame = tk.LabelFrame(entry_frame, text="")
    furtherproc_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=30, pady=15)

    furtherproc_var = tk.IntVar(value=True)

    furtherproc_entry = tk.Checkbutton(
        furtherproc_frame,
        text="  I want to cyclize or further modify my peptide",
        variable=furtherproc_var,
        onvalue=True,
        offvalue=False,
    )
    furtherproc_entry.grid(row=0, column=0, sticky="W")

    for widget in furtherproc_frame.winfo_children():
        widget.grid_configure(padx=20, pady=5)

    button = tk.Button(
        entry_frame,
        text="CREATE PEPTIDE",
        bg="#0052cc",
        fg="#ffffff",
        command=userinput1_click_enter_button,
    )
    button.pack(padx=20, pady=(15, 30))

    aminoacid_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aminoacid_frame.grid(row=0, column=1, padx=50, pady=100, sticky="news")
    df_aminoacids_display = pd.DataFrame(
        data={
            "Title": [
                "Natural",
                "D and N-Me",
                "Non-natural",
                "Peptoid",
                "Own Amino Acids",
            ],
            "Mols": [[], [], [], [], []],
        }
    )

    aalibrarytitle_label = tk.Label(
        aminoacid_frame, text="AVAILABLE AMINO ACIDS", bg="#d9d9d9"
    )
    aalibrarytitle_label.pack(side=tk.TOP, fill=tk.BOTH, padx=30, pady=(30, 0))

    def click_in_entry(event, entry_field, default_text):
        if entry_field.get() == default_text:
            entry_field.delete(0, tk.END)
            entry_field.config(fg="black")

    def click_out_entry(event, entry_field, default_text):
        if not entry_field.get():
            entry_field.insert(0, default_text)
            entry_field.config(fg="#949494")

    aalibrary_frame = tk.LabelFrame(aminoacid_frame, text="")
    aalibrary_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=30, pady=15)

    default_aa_frame = tk.Frame(aalibrary_frame)
    default_aa_frame.pack(side=tk.TOP, padx=10)

    option_lst = df_aminoacids_display.loc[:3, "Title"].tolist()
    for idx, option in enumerate(option_lst):
        if idx == 0 or idx == 3:
            w = 7
        else:
            w = 10
        aadefault_button = tk.Button(
            default_aa_frame,
            text=option,
            width=w,
            bg="#949494",
            fg="#ffffff",
            command=lambda t=option: userinput1_click_aabutton(t),
        )
        aadefault_button.grid(row=1, column=idx, padx=8, pady=15, sticky="news")

    aalibrarytitle_label = tk.Label(
        aminoacid_frame, text="ADD OWN AMINO ACIDS", bg="#d9d9d9"
    )
    aalibrarytitle_label.pack(side=tk.TOP, fill=tk.BOTH, padx=30, pady=(10, 0))

    ownaa_frame = tk.LabelFrame(aminoacid_frame, text="")
    ownaa_frame.pack(fill=tk.BOTH, expand=True, padx=30, pady=15)

    addaa_frame = tk.Frame(ownaa_frame)
    addaa_frame.pack(side=tk.TOP, pady=10)
    addaa_entrycode_entry = tk.Entry(addaa_frame, width=15, fg="#949494")
    addaa_entrycode_entry.insert(0, "Entry Code")
    addaa_entrycode_entry.bind(
        "<FocusIn>",
        lambda event: click_in_entry(event, addaa_entrycode_entry, "Entry Code"),
    )
    addaa_entrycode_entry.bind(
        "<FocusOut>",
        lambda event: click_out_entry(event, addaa_entrycode_entry, "Entry Code"),
    )
    addaa_entrycode_entry.grid(row=1, column=0, padx=5, pady=(10, 5))
    addaa_smiles_entry = tk.Entry(addaa_frame, width=25, fg="#949494")
    addaa_smiles_entry.insert(0, "Smiles Code")
    addaa_smiles_entry.bind(
        "<FocusIn>",
        lambda event: click_in_entry(event, addaa_smiles_entry, "Smiles Code"),
    )
    addaa_smiles_entry.bind(
        "<FocusOut>",
        lambda event: click_out_entry(event, addaa_smiles_entry, "Smiles Code"),
    )
    addaa_smiles_entry.grid(row=1, column=1, padx=5, pady=(10, 5))
    addaa_button = tk.Button(
        addaa_frame,
        text=" + ",
        bg="#949494",
        fg="#ffffff",
        command=lambda: click_addaa(),
    )
    addaa_button.grid(row=1, column=2, padx=5, pady=(10, 5))

    addaa_buttons_frame = tk.Frame(ownaa_frame)
    addaa_buttons_frame.pack(side=tk.TOP, pady=(0, 15))
    ownaa_button = tk.Button(
        addaa_buttons_frame,
        text="Import List",
        width=10,
        command=lambda: userinput1_click_importaabutton(),
    )
    ownaa_button.grid(row=0, column=0, padx=9, pady=5, sticky="news")
    ownaa_button2 = tk.Button(
        addaa_buttons_frame,
        text="Export Added",
        width=10,
        command=lambda: userinput1_click_exportaabutton(),
    )
    ownaa_button2.grid(row=0, column=1, padx=9, pady=5, sticky="news")
    ownaa_button3 = tk.Button(
        addaa_buttons_frame,
        text="Overview",
        width=8,
        command=lambda t="Own Amino Acids": userinput1_click_aabutton(t),
    )
    ownaa_button3.grid(row=0, column=2, padx=9, pady=5, sticky="news")
    info_button2 = tk.Button(
        addaa_buttons_frame,
        image=info_button_png,
        command=help_button1_aminoacids_click,
        bd=0,
        highlightthickness=0,
    )
    info_button2.grid(row=0, column=3, padx=(2, 13), pady=5)

    aalibrarytitle_label = tk.Label(
        aminoacid_frame, text="MULTIPLE PEPTIDE CALCULATION", bg="#d9d9d9"
    )
    aalibrarytitle_label.pack(side=tk.TOP, fill=tk.BOTH, padx=30, pady=(10, 0))

    fastcalc_frame = tk.LabelFrame(aminoacid_frame, text="")
    fastcalc_frame.pack(fill=tk.BOTH, expand=True, padx=30, pady=(15, 30))

    fastcalc_buttons_frame = tk.Frame(fastcalc_frame)
    fastcalc_buttons_frame.pack(side=tk.TOP)

    startfastcalc_button = tk.Button(
        fastcalc_buttons_frame,
        text="START CALCULATION",
        bg="#0052cc",
        fg="#ffffff",
        width=20,
        command=lambda: click_startcalculation_button(),
    )  #'#0052cc''#949494'
    startfastcalc_button.grid(row=0, column=0, padx=10, pady=20, sticky="news")
    info_button1 = tk.Button(
        fastcalc_buttons_frame,
        image=info_button_png,
        command=help_button1_multiplepepcalc_click,
        bd=0,
        highlightthickness=0,
    )  # , bg = white)
    info_button1.grid(row=0, column=1, padx=10, pady=15, sticky="news")

    icon_frame = tk.Frame(window_frame)
    icon_frame.grid(row=2, columnspan=2, pady=10)
    teriyaki_icon = tk.Label(icon_frame, image=Teriyaki)
    teriyaki_icon.grid(row=0, column=0)
    help_png = Image.open(path_to_helpbutton)
    resize_help_png = help_png.resize((25, 25))
    help_button_png = ImageTk.PhotoImage(resize_help_png)
    help_button = tk.Button(
        icon_frame,
        image=help_button_png,
        command=help_button1_click,
        bd=0,
        highlightthickness=0,
    )
    help_button.grid(row=0, column=1)

    window_size_centered(window, 1422, 800)
    window.iconphoto(False, Teriyaki)

    window.bind("<Return>", userinput1_click_enter_button)
    window.bind("<Escape>", click_close_button)
    window.bind("<question>", help_button1_click)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.mainloop()

    return (
        sequence,
        cterm,
        nterm,
        cycl,
        window_status,
        default_pept,
        df_aminoacids,
        calc_type,
        calc_table,
        calc_table_file_path,
    )


## USER INPUT 2 ##
## Cyclization and Modification choice ##


def userinput2(
    cyclization_choices, aamod_choices, aa_list, peptide_memory, window_status, path_to_teriyaki, path_to_teriyaki_easteregg
):

    def userinput2_click_cycle_button():
        nonlocal cyclization_finalchoice, cycl, aamod
        cyclization_finalchoice = cyclization_selected_choice.get()
        cycl = True
        aamod = False
        main_window2.destroy()

    def userinput2_click_modify_button():
        nonlocal aamod_finalchoice, aamod, cycl
        aamod_finalchoice = aamod_selected_choice.get()
        aamod = True
        cycl = False
        main_window2.destroy()

    def userinput2_click_done_button(event=None):
        nonlocal window_status
        window_status = 3
        main_window2.destroy()

    def userinput2_click_startagain_button(event=None):
        nonlocal window_status
        window_status = 1
        main_window2.destroy()

    def userinput2_click_back_button(event=None):
        nonlocal back, peptide_memory, window_status
        back = True
        window_status = 2
        peptide_memory = peptide_memory[:-1]
        main_window2.destroy()

    def click_close_button(event=None):
        nonlocal window_status
        done_input = messagebox.askquestion(
            " ",
            "Are you sure you want to close the program? \n \nATTENTION: Self-added amino acids will be lost! \nSave custom amino acids before closing the program.",
        )
        if done_input == "yes":
            window_status = 4
            main_window2.destroy()

    def click_getsmiles_button():
        nonlocal peptide_memory
        smiles = Chem.MolToSmiles(peptide_memory[-1].peptide)
        pyperclip.copy(smiles)
        messagebox.showinfo("Success", "Value copied to clipboard")

    def click_savejpg_button():
        nonlocal molecule
        save_jpg_path = filedialog.asksaveasfilename(
            defaultextension=".jpg",
            filetypes=[("JPEG files", "*.jpg"), ("All files", "*.*")],
        )
        if save_jpg_path:
            watermark = "Created with PICKAPEP"
            current_peptide_image_save = Draw.MolToImage(molecule, size=(700, 412))
            current_peptide_draw = ImageDraw.Draw(current_peptide_image_save)
            text_font = ImageFont.truetype("arial.ttf", 12)
            text_width, text_height = current_peptide_draw.textsize(
                watermark, font=text_font
            )
            text_margin = 10
            text_position = (
                current_peptide_image_save.width - text_width - text_margin,
                current_peptide_image_save.height - text_height - text_margin,
            )
            current_peptide_draw.text(
                text_position, watermark, fill=(179, 179, 179), font=text_font
            )
            current_peptide_image_save.save(save_jpg_path, dpi=(300, 300))
            messagebox.showinfo("Save", " Image saved successfully")

    def click_getsummary_properties():
        nonlocal df_summary_properties
        save_file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=(
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ),
        )
        if save_file_path:
            if save_file_path.endswith(".csv"):
                df_summary_properties.to_csv(save_file_path, index=False)
            elif save_file_path.endswith(".xlsx"):
                df_summary_properties.to_excel(save_file_path, index=False)
            messagebox.showinfo("Save", "File saved successfully")

    def click_get3Dpymol_button():
        nonlocal peptide_memory
        generate_threedim_pymol(peptide_memory)

    def help_button2_click(event=None):
        subwindow_help2 = tk.Toplevel(main_window2)
        subwindow_help2.title("PICKAPEP SUPPORT")

        help2_display = tk.Frame(subwindow_help2, bg="#ffffff")
        help2_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help2.iconphoto(False, Teriyaki)

        img_Help2 = Image.open(path_to_HS_Userinput2)
        img_Help2_width, img_Help2_height = img_Help2.size
        ratio_size = img_Help2_height / img_Help2_width
        target_width = 800
        img_Help2 = img_Help2.resize((target_width, int(target_width * ratio_size)))
        img_Help2_tk = ImageTk.PhotoImage(img_Help2)
        Help2_display = tk.Label(help2_display, image=img_Help2_tk, bg="#ffffff")
        Help2_display.pack(padx=50, pady=50)

        subwindow_help2.mainloop()

    def help_button2_cyclization_click():

        subwindow_help2_cycl = tk.Toplevel(main_window2)
        subwindow_help2_cycl.title("PICKAPEP SUPPORT")

        help2_cycl_display = tk.Frame(subwindow_help2_cycl, bg="#ffffff")
        help2_cycl_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help2_cycl.iconphoto(False, Teriyaki)

        img_Help2_cycl = Image.open(path_to_HS_Userinput2_cycl)
        img_Help2_cycl_width, img_Help2_cycl_height = img_Help2_cycl.size
        ratio_size = img_Help2_cycl_height / img_Help2_cycl_width
        target_width = 800
        img_Help2_cycl = img_Help2_cycl.resize(
            (target_width, int(target_width * ratio_size))
        )
        img_Help2_cycl_tk = ImageTk.PhotoImage(img_Help2_cycl)
        Help2_cycl_display = tk.Label(
            help2_cycl_display, image=img_Help2_cycl_tk, bg="#ffffff"
        )
        Help2_cycl_display.pack(padx=50, pady=50)

        subwindow_help2_cycl.mainloop()

    def help_button2_modification_click():

        subwindow_help2_mod = tk.Toplevel(main_window2)
        subwindow_help2_mod.title("PICKAPEP SUPPORT")

        help2_mod_display = tk.Frame(subwindow_help2_mod, bg="#ffffff")
        help2_mod_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help2_mod.iconphoto(False, Teriyaki)

        img_Help2_mod = Image.open(path_to_HS_Userinput2_mod)
        img_Help2_mod_width, img_Help2_mod_height = img_Help2_mod.size
        ratio_size = img_Help2_mod_height / img_Help2_mod_width
        target_width = 800
        img_Help2_mod = img_Help2_mod.resize(
            (target_width, int(target_width * ratio_size))
        )
        img_Help2_mod_tk = ImageTk.PhotoImage(img_Help2_mod)
        Help2_mod_display = tk.Label(
            help2_mod_display, image=img_Help2_mod_tk, bg="#ffffff"
        )
        Help2_mod_display.pack(padx=50, pady=50)

        subwindow_help2_mod.mainloop()

    cycl = False
    cyclization_finalchoice = ""
    aamod = False
    aamod_finalchoice = ""
    back = False
    curr_pep = peptide_memory[-1]
    molecule = curr_pep.peptide
    df_summary_properties = create_peptide_summary(peptide_memory)

    main_window2 = tk.Tk()
    main_window2.title("PICKAPEP")

    sequence_easteregg = "".join(
        list(peptide_memory[-1].df_sequence["AminoAcid"][1:-1])
    )
    if sequence_easteregg == "TERIYAKI":
        path_to_teriyaki = path_to_teriyaki_easteregg
        txt_done = "EGG-xellent"
    else:
        txt_done = "Done"

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    main_window2_frame = tk.Frame(main_window2)
    main_window2_frame.pack()

    main_window2_frame.columnconfigure(0, minsize=550)
    main_window2_frame.columnconfigure(1, minsize=550)

    current_structure_label_text = "YOUR CURRENT PEPTIDE"
    overview_label_text = "MODIFICATION OVERVIEW"
    detail_label_text = "PEPTIDE INFORMATION"
    background_colour = "#d9d9d9"

    win1_current_peptide_frame = tk.LabelFrame(
        main_window2_frame, text="", bg=background_colour
    )
    win1_current_peptide_frame.grid(
        row=0, column=0, sticky="news", padx=(20, 10), pady=(25, 10)
    )

    current_structure_label = tk.Label(
        win1_current_peptide_frame,
        text=current_structure_label_text,
        bg=background_colour,
    )
    current_structure_label.pack(padx=10, pady=(15, 0))

    current_peptide_display = tk.LabelFrame(win1_current_peptide_frame, text="")
    current_peptide_display.pack(padx=15, pady=(10, 0))
    rdCoordGen.AddCoords(molecule)
    current_peptide_image = Draw.MolToImage(molecule, size=(465, 250))
    current_peptide_photo = ImageTk.PhotoImage(current_peptide_image)
    current_peptide_photo_display = tk.Label(
        current_peptide_display, image=current_peptide_photo
    )
    current_peptide_photo_display.pack()

    getinfo_frame_sub = tk.Frame(win1_current_peptide_frame, bg=background_colour)
    getinfo_frame_sub.pack(padx=15, pady=0)

    getsmiles_button = tk.Button(
        getinfo_frame_sub,
        text="Get Smiles",
        bg="#949494",
        fg="#ffffff",
        width=10,
        command=lambda: click_getsmiles_button(),
    )
    getsmiles_button.grid(row=0, column=0, padx=10, pady=(10, 15))
    savejpg_button = tk.Button(
        getinfo_frame_sub,
        text="Save Image",
        bg="#949494",
        fg="#ffffff",
        width=10,
        command=lambda: click_savejpg_button(),
    )
    savejpg_button.grid(row=0, column=1, padx=10, pady=(10, 15))
    getinfo_button = tk.Button(
        getinfo_frame_sub,
        text="Peptide Summary",
        bg="#949494",
        fg="#ffffff",
        width=15,
        command=lambda: click_getsummary_properties(),
    )
    getinfo_button.grid(row=0, column=2, padx=10, pady=(10, 15))
    get3Dpymol_button = tk.Button(
        getinfo_frame_sub,
        text="PyMOL 3D View",
        bg="#949494",
        fg="#ffffff",
        width=15,
        command=lambda: click_get3Dpymol_button(),
    )
    get3Dpymol_button.grid(row=0, column=3, padx=10, pady=(10, 15))

    win2_entry_frame = tk.LabelFrame(main_window2_frame, text="", bg=background_colour)
    win2_entry_frame.grid(row=0, column=1, sticky="news", padx=(10, 20), pady=(25, 10))

    win2_entry_frame.columnconfigure(0, weight=1)
    win2_entry_frame.columnconfigure(1, weight=1)

    info_png = Image.open(path_to_infobutton)
    resize_info_png = info_png.resize((20, 20))
    info_button_png = ImageTk.PhotoImage(resize_info_png)

    options_label = tk.Label(
        win2_entry_frame, text="MODIFICATION OPTIONS", bg=background_colour
    )
    options_label.grid(
        row=0, column=0, columnspan=2, padx=10, pady=(15, 0), sticky="news"
    )

    cyclization_label_frame = tk.Frame(win2_entry_frame, bg=background_colour)
    cyclization_label_frame.grid(row=1, column=0, padx=10, pady=(15, 0), sticky="news")

    cyclization_label = tk.Label(
        cyclization_label_frame, text="Cyclization Reactions", bg=background_colour
    )
    cyclization_label.pack(side=tk.LEFT, padx=(55, 0), pady=(0, 10))

    info_button1 = tk.Button(
        cyclization_label_frame,
        image=info_button_png,
        command=help_button2_cyclization_click,
        bd=0,
        highlightthickness=0,
        bg=background_colour,
    )
    info_button1.pack(side=tk.LEFT, padx=5, pady=(0, 10))

    cyclization_frame = tk.LabelFrame(win2_entry_frame, text="")
    cyclization_frame.grid(row=2, column=0, padx=(30, 15), pady=5, sticky="news")

    cyclization_selected_choice = tk.StringVar(value=cyclization_choices[0])

    for idx, choice in enumerate(cyclization_choices):
        radiobutton_cyclization = tk.Radiobutton(
            cyclization_frame,
            text=choice,
            variable=cyclization_selected_choice,
            value=choice,
        )
        if idx == 0:
            radiobutton_cyclization.grid(
                row=idx, column=0, padx=10, pady=(10, 3), sticky="w"
            )
        elif idx != len(cyclization_choices) - 1:
            radiobutton_cyclization.grid(row=idx, column=0, padx=10, pady=3, sticky="w")
        else:
            radiobutton_cyclization.grid(
                row=idx, column=0, padx=10, pady=(3, 10), sticky="w"
            )

    button_cyclization = tk.Button(
        win2_entry_frame,
        text="Start Cyclization",
        bg="#0052cc",
        fg="#ffffff",
        command=lambda: userinput2_click_cycle_button(),
    )
    button_cyclization.grid(row=3, column=0, padx=(30, 15), pady=(15, 5))

    aamod_label_frame = tk.Frame(win2_entry_frame, bg=background_colour)
    aamod_label_frame.grid(row=1, column=1, padx=10, pady=(15, 0), sticky="news")

    aamod_label = tk.Label(
        aamod_label_frame, text="Amino Acid Modifications", bg=background_colour
    )
    aamod_label.pack(side=tk.LEFT, padx=(35, 0), pady=(0, 10))

    info_button2 = tk.Button(
        aamod_label_frame,
        image=info_button_png,
        command=help_button2_modification_click,
        bd=0,
        highlightthickness=0,
        bg=background_colour,
    )
    info_button2.pack(side=tk.LEFT, padx=5, pady=(0, 10))

    aamod_frame = tk.LabelFrame(win2_entry_frame, text="")
    aamod_frame.grid(row=2, column=1, padx=(15, 30), pady=5, sticky="news")

    aamod_selected_choice = tk.StringVar(value=aamod_choices[0])

    for idx, choice in enumerate(aamod_choices):
        radiobutton_aamod = tk.Radiobutton(
            aamod_frame, text=choice, variable=aamod_selected_choice, value=choice
        )
        if idx == 0:
            radiobutton_aamod.grid(row=idx, column=0, padx=10, pady=(10, 3), sticky="w")
        elif idx != len(aamod_choices) - 1:
            radiobutton_aamod.grid(row=idx, column=0, padx=10, pady=3, sticky="w")
        else:
            radiobutton_aamod.grid(row=idx, column=0, padx=10, pady=(3, 10), sticky="w")

    button_aamod = tk.Button(
        win2_entry_frame,
        text="Start Modification",
        bg="#0052cc",
        fg="#ffffff",
        command=lambda: userinput2_click_modify_button(),
    )
    button_aamod.grid(row=3, column=1, padx=(15, 30), pady=(15, 5))

    win3_overview_frame = tk.LabelFrame(
        main_window2_frame, text="", bg=background_colour
    )
    win3_overview_frame.grid(
        row=1, column=0, sticky="news", padx=(20, 10), pady=(10, 20)
    )

    overview_label = tk.Label(
        win3_overview_frame, text=overview_label_text, bg=background_colour
    )
    overview_label.pack(padx=10, pady=(15, 0))

    overview_frame_sub = tk.LabelFrame(win3_overview_frame, text="", bg="#FFFFFF")
    overview_frame_sub.columnconfigure(0, minsize=130)
    overview_frame_sub.columnconfigure(1, minsize=300)
    overview_frame_sub.pack(pady=(10, 20))

    curr_sequence_frame = tk.Frame(overview_frame_sub, bg="#FFFFFF")
    curr_sequence_frame.grid(
        row=0, column=0, rowspan=2, sticky="news", padx=10, pady=(8, 18)
    )

    curr_sequence_label = tk.Label(curr_sequence_frame, text="Sequence", bg="#FFFFFF")
    curr_sequence_label.pack(padx=10, pady=(0, 5))
    curr_sequence_listbox = tk.Listbox(curr_sequence_frame, height=9, width=20)

    curr_sequence_listbox.pack(side="left", fill="both")
    curr_sequence_scrollbar = tk.Scrollbar(curr_sequence_frame, troughcolor="#d9d9d9")
    curr_sequence_scrollbar.pack(side="right", fill="both")

    for aa in aa_list:
        curr_sequence_listbox.insert("end", aa)
    for i in range(0, len(aa_list), 2):
        curr_sequence_listbox.itemconfig(i, bg="#f0f0f0")

    curr_sequence_listbox.config(yscrollcommand=curr_sequence_scrollbar.set)
    curr_sequence_scrollbar.config(command=curr_sequence_listbox.yview)

    curr_cycl_frame = tk.Frame(overview_frame_sub, bg="#FFFFFF")
    curr_cycl_frame.grid(row=0, column=1, sticky="news", padx=10, pady=(8, 4))

    curr_cycl_label = tk.Label(curr_cycl_frame, text="Cyclization", bg="#FFFFFF")
    curr_cycl_label.pack(padx=10, pady=(0, 5))
    curr_cycl_listbox = tk.Listbox(curr_cycl_frame, height=4, width=45)
    curr_cycl_listbox.pack(side="left", fill="both")
    curr_cycl_scrollbar = tk.Scrollbar(curr_cycl_frame, troughcolor="#d9d9d9")
    curr_cycl_scrollbar.pack(side="right", fill="both")

    cycl_info_list, mod_info_list = get_info_lists(peptide_memory)

    if len(cycl_info_list) < 4:
        elements_needed = 4 - len(cycl_info_list)
        for i in range(elements_needed):
            cycl_info_list.append(" ")

    for cycl_info in cycl_info_list:
        curr_cycl_listbox.insert("end", cycl_info)
    for i in range(0, len(cycl_info_list), 2):
        curr_cycl_listbox.itemconfig(i, bg="#f0f0f0")

    curr_cycl_listbox.config(yscrollcommand=curr_cycl_scrollbar.set)
    curr_cycl_scrollbar.config(command=curr_cycl_listbox.yview)

    curr_mod_frame = tk.Frame(overview_frame_sub, bg="#FFFFFF")
    curr_mod_frame.grid(row=1, column=1, sticky="news", padx=10, pady=(4, 18))

    curr_mod_label = tk.Label(curr_mod_frame, text="Modification", bg="#FFFFFF")
    curr_mod_label.pack(padx=10, pady=(0, 5))
    curr_mod_listbox = tk.Listbox(curr_mod_frame, height=4, width=45)
    curr_mod_listbox.pack(side="left", fill="both")
    curr_mod_scrollbar = tk.Scrollbar(curr_mod_frame, troughcolor="#d9d9d9")
    curr_mod_scrollbar.pack(side="right", fill="both")

    if len(mod_info_list) < 4:
        elements_needed = 4 - len(mod_info_list)
        for i in range(elements_needed):
            mod_info_list.append(" ")

    for mod_info in mod_info_list:
        curr_mod_listbox.insert("end", mod_info)
    for i in range(0, len(mod_info_list), 2):
        curr_mod_listbox.itemconfig(i, bg="#f0f0f0")

    curr_mod_listbox.config(yscrollcommand=curr_mod_scrollbar.set)
    curr_mod_scrollbar.config(command=curr_mod_listbox.yview)

    win4_detail_frame = tk.LabelFrame(main_window2_frame, text="", bg=background_colour)
    win4_detail_frame.grid(row=1, column=1, sticky="news", padx=(10, 20), pady=(10, 20))

    detail_label = tk.Label(
        win4_detail_frame, text=detail_label_text, bg=background_colour
    )
    detail_label.pack(padx=10, pady=(15, 0))

    detail_frame_sub = tk.LabelFrame(win4_detail_frame, text="", bg="#FFFFFF")
    detail_frame_sub.pack(padx=15, pady=(10, 15))

    detail_frame_sub2 = tk.Frame(detail_frame_sub, bg="#FFFFFF")
    detail_frame_sub2.columnconfigure(0, minsize=240)
    detail_frame_sub2.columnconfigure(1, minsize=220)
    detail_frame_sub2.pack(pady=10)
    df_details = curr_pep.pepinfo

    for i, row in df_details.iterrows():
        detail_text = f"{row['Info']}:   {row['Value']}"
        if i < 7:
            width_label = 40
        else:
            width_label = 30
        detail_text_label = tk.Text(
            detail_frame_sub2,
            bg="#FFFFFF",
            width=width_label,
            height=1,
            bd=0,
            highlightthickness=0,
        )
        detail_text_label.configure(font=("Segoe UI", 9))
        detail_text_label.insert(tk.END, detail_text)
        pos_colon = detail_text.find(":")
        detail_text_label.tag_add("italic", "1.0", f"1.{pos_colon}")
        detail_text_label.tag_configure("italic", font=("Segoe UI", 9, "normal"))
        if row["Info"] == "Chemical Formula":
            for letter, char in enumerate(detail_text):
                if char.isdigit():
                    detail_text_label.tag_add(
                        "subscript", f"1.{letter}", f"1.{letter+1}"
                    )
            detail_text_label.tag_configure(
                "subscript", offset=-1, font=("Segoe UI", 7, "normal")
            )
        elif detail_text[5] == "+":
            detail_text_label.tag_add("superscript", f"1.{5}", f"1.{6}")
            detail_text_label.tag_configure(
                "superscript", offset=4, font=("Segoe UI", 6, "normal")
            )
        elif detail_text[7] == "+":
            detail_text_label.tag_add("superscript", f"1.{6}", f"1.{8}")
            detail_text_label.tag_configure(
                "superscript", offset=4, font=("Segoe UI", 6, "normal")
            )

        if i < 7:
            detail_text_label.grid(row=i, column=0, padx=(15, 0), pady=6, sticky="w")
        else:
            detail_text_label.grid(row=i - 7, column=1, padx=(0, 0), pady=6, sticky="w")
        detail_text_label.config(state=tk.DISABLED)

    done_frame = tk.Frame(main_window2_frame)
    done_frame.grid(row=2, column=0, columnspan=2)
    startagain_button = tk.Button(
        done_frame,
        text="Start Again",
        bg="#d9d9d9",
        width=10,
        command=lambda: userinput2_click_startagain_button(),
    )
    startagain_button.grid(row=0, column=0)
    if len(peptide_memory) > 1:
        back_button = tk.Button(
            done_frame,
            text="Undo",
            bg="#d9d9d9",
            width=10,
            command=lambda: userinput2_click_back_button(),
        )
        back_button.grid(row=0, column=2)
        main_window2.bind("<BackSpace>", userinput2_click_back_button)
    done_button = tk.Button(
        done_frame,
        text=txt_done,
        bg="#d9d9d9",
        width=10,
        command=lambda: userinput2_click_done_button(),
    )
    done_button.grid(row=0, column=3)

    for widget in done_frame.winfo_children():
        widget.grid_configure(padx=5, pady=(5, 0))

    icon_frame = tk.Frame(main_window2_frame)
    icon_frame.grid(row=3, columnspan=2, pady=5)
    teriyaki_icon = tk.Label(icon_frame, image=Teriyaki)
    teriyaki_icon.grid(row=0, column=0)
    help_png = Image.open(path_to_helpbutton)
    resize_help_png = help_png.resize((25, 25))
    help_button_png = ImageTk.PhotoImage(resize_help_png)
    help_button = tk.Button(
        icon_frame,
        image=help_button_png,
        command=help_button2_click,
        bd=0,
        highlightthickness=0,
    )
    help_button.grid(row=0, column=1)

    window_size_centered(main_window2, 1422, 800)
    main_window2.iconphoto(False, Teriyaki)

    main_window2.protocol("WM_DELETE_WINDOW", click_close_button)

    main_window2.bind("<Return>", userinput2_click_done_button)
    main_window2.bind("<Escape>", click_close_button)
    main_window2.bind("<Delete>", userinput2_click_startagain_button)

    main_window2.bind("<question>", help_button2_click)

    main_window2.mainloop()

    return (
        cycl,
        cyclization_finalchoice,
        aamod,
        aamod_finalchoice,
        peptide_memory,
        back,
        window_status,
    )


## USER INPUT 2 SUBWINDOW Cyclization##


def userinput2_cycl(
    cyclization_finalchoice,
    aa_options_input,
    additional_needed,
    additional_options,
    additional_info,
    path_to_teriyaki,
):

    def click_close_button(event=None):
        nonlocal skip_cycl
        skip_cycl = True
        window.destroy()

    def userinput2_cycl_click_enter_button(event=None):
        nonlocal aa_choice_output, additional_output

        aa_choice_output = [var.get() for var in variable_entry_list]

        if additional_needed:
            additional_output = additional_combobox.get()

        if all(
            any(element in option for option in aa_options_input)
            for element in aa_choice_output
        ) and all(element != "" for element in aa_choice_output):
            if additional_needed:
                if additional_output in additional_options and additional_output != "":
                    window.destroy()
                else:
                    tk.messagebox.showwarning(
                        message="Please choose a valid " + additional_info
                    )
            else:
                window.destroy()
        else:
            tk.messagebox.showwarning(message="Please choose valid amino acids")

    aa_choice_output = []
    additional_output = 0
    skip_cycl = False

    window = tk.Tk()
    window.title("Cyclization")

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    window_frame = tk.Frame(window)
    window_frame.pack(padx=20, pady=(30, 30))

    aa_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aa_frame.pack(pady=(60, 0))

    aa_label1 = tk.Label(aa_frame, text=cyclization_finalchoice.upper(), bg="#d9d9d9")
    aa_label1.pack(padx=20, pady=(30, 0))

    cyclaa_entry_frame = tk.LabelFrame(aa_frame, text="", bg="#f0f0f0")
    cyclaa_entry_frame.pack(padx=70, pady=20)

    variable_entry_list = []

    aachoice_labeltxt = tk.Label(
        cyclaa_entry_frame, text="Choose Amino Acids", bg="#f0f0f0"
    )
    aachoice_labeltxt.grid(
        row=0, columnspan=len(aa_options_input), padx=10, pady=(15, 0)
    )

    for i, option in enumerate(aa_options_input):
        aachoice_var = tk.StringVar(cyclaa_entry_frame)
        aachoice_var.set(option[0])
        if i == 0:
            x_value = (20, 10)
        elif i == len(aa_options_input) - 1:
            x_value = (10, 20)
        else:
            x_value = 10
        combobox = ttk.Combobox(
            cyclaa_entry_frame, textvariable=aachoice_var, values=option, width=15
        )
        combobox.grid(row=1, column=i, padx=x_value, pady=(10, 20))
        variable_entry_list.append(aachoice_var)

    if additional_needed:
        additional_entry_frame = tk.LabelFrame(aa_frame, text="")
        additional_entry_frame.pack(side=tk.TOP, fill=tk.BOTH, padx=70, pady=(0, 20))
        additional_label = tk.Label(
            additional_entry_frame, text="Choose the Additional Component"
        )
        additional_label.pack(padx=10, pady=(15, 0))

        additional_combobox = ttk.Combobox(
            additional_entry_frame, values=additional_options, width=30
        )
        additional_combobox.pack(padx=30, pady=(10, 20))

    button_frame = tk.Frame(aa_frame, bg="#d9d9d9")
    button_frame.pack()

    button_back = tk.Button(
        button_frame, text="Back", bg="#f0f0f0", command=click_close_button, width=8
    )
    button_back.grid(row=0, column=0, padx=5, pady=(15, 30))

    button = tk.Button(
        button_frame,
        text="START",
        bg="#0052cc",
        fg="#ffffff",
        command=userinput2_cycl_click_enter_button,
        width=8,
    )
    button.grid(row=0, column=1, padx=5, pady=(15, 30))

    window_size_centered(window, 1067, 600)
    window.iconphoto(False, Teriyaki)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.bind("<Return>", userinput2_cycl_click_enter_button)
    window.bind("<Escape>", click_close_button)

    window.mainloop()
    return aa_choice_output, additional_output, skip_cycl


## USER INPUT 2 SUBWINDOW Fatty acid modification##


def userinput2_aamod(aa_options, aamod_finalchoice, path_to_teriyaki):

    def click_close_button(event=None):
        nonlocal aa_choice
        aa_choice = ""
        window.destroy()

    def userinput2_aamod_click_enter_button(event=None):
        nonlocal aa_choice
        aa_choice = aa_combobox.get()
        if aa_choice in aa_options and aa_choice != "":
            window.destroy()
        else:
            tk.messagebox.showwarning(message="Please choose valid option")

    aa_choice = ""

    window = tk.Tk()
    window.title("Amino Acid Modification")

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    window_frame = tk.Frame(window)
    window_frame.pack(padx=20, pady=(30, 30))

    aa_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aa_frame.pack(padx=20, pady=(60, 0))

    aa_label1 = tk.Label(aa_frame, text=aamod_finalchoice.upper(), bg="#d9d9d9")
    aa_label1.pack(padx=20, pady=(30, 0))

    aamodaa_entry_frame = tk.LabelFrame(aa_frame, text="", bg="#f0f0f0")
    aamodaa_entry_frame.pack(padx=70, pady=20)

    aachoice_labeltxt = tk.Label(
        aamodaa_entry_frame, text="Choose Amino Acid", bg="#f0f0f0"
    )
    aachoice_labeltxt.pack(padx=10, pady=(15, 0))

    aa_combobox = ttk.Combobox(aamodaa_entry_frame, values=aa_options)
    aa_combobox.set(aa_options[0])
    aa_combobox.pack(padx=20, pady=(10, 20))

    button_frame = tk.Frame(aa_frame, bg="#d9d9d9")
    button_frame.pack()

    button_back = tk.Button(
        button_frame, text="Back", bg="#f0f0f0", command=click_close_button, width=8
    )
    button_back.grid(row=0, column=0, padx=5, pady=(15, 30))

    button = tk.Button(
        button_frame,
        text="START",
        bg="#0052cc",
        fg="#ffffff",
        command=userinput2_aamod_click_enter_button,
        width=8,
    )
    button.grid(row=0, column=1, padx=5, pady=(15, 30))

    window_size_centered(window, 1067, 600)
    window.iconphoto(False, Teriyaki)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.bind("<Return>", userinput2_aamod_click_enter_button)
    window.bind("<Escape>", click_close_button)

    window.mainloop()

    return aa_choice


def userinput2_aamod_option(aa_options, choices, aamod_finalchoice, path_to_teriyaki):

    def click_close_button(event=None):
        nonlocal aa_choice
        aa_choice = ""
        window.destroy()

    def userinput2_aamod_option_click_enter_button(event=None):
        nonlocal aa_choice, finalchoice
        finalchoice = selected_choice.get()
        aa_choice = aa_combobox.get()
        if aa_choice in aa_options and aa_choice != "":
            window.destroy()
        else:
            tk.messagebox.showwarning(message="Please choose valid option")

    aa_choice = ""
    finalchoice = []

    window = tk.Tk()
    window.title("Amino Acid Modification")

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    window_frame = tk.Frame(window)
    window_frame.pack(padx=20, pady=(30, 30))

    aa_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aa_frame.pack(padx=20, pady=(30, 0))

    aa_label1 = tk.Label(aa_frame, text=aamod_finalchoice.upper(), bg="#d9d9d9")
    aa_label1.pack(padx=20, pady=(30, 0))

    aamodaa_entry_frame = tk.LabelFrame(aa_frame, text="", bg="#f0f0f0")
    aamodaa_entry_frame.pack(padx=70, pady=20)

    aachoice_labeltxt = tk.Label(
        aamodaa_entry_frame, text="Choose Amino Acid", bg="#f0f0f0"
    )
    aachoice_labeltxt.pack(padx=10, pady=(15, 0))

    aa_combobox = ttk.Combobox(aamodaa_entry_frame, values=aa_options)
    aa_combobox.set(aa_options[0])
    aa_combobox.pack(padx=20, pady=(10, 20))

    choices_frame = tk.LabelFrame(aa_frame, text="")
    choices_frame.pack(padx=60, pady=(10, 20))

    choices_frame_sub = tk.Frame(choices_frame)
    choices_frame_sub.pack(padx=40, pady=20)

    selected_choice = tk.StringVar(value=choices[0])

    for idx, choice in enumerate(choices):
        radiobutton = tk.Radiobutton(
            choices_frame_sub, text=choice, variable=selected_choice, value=choice
        )
        if idx < 6:
            radiobutton.grid(row=idx, column=0, sticky="W", padx=(20, 30))
        else:
            radiobutton.grid(row=idx - 6, column=1, sticky="W", padx=30)

    button_frame = tk.Frame(aa_frame, bg="#d9d9d9")
    button_frame.pack()

    button_back = tk.Button(
        button_frame, text="Back", bg="#f0f0f0", command=click_close_button, width=8
    )
    button_back.grid(row=0, column=0, padx=5, pady=(5, 30))

    button = tk.Button(
        button_frame,
        text="START",
        bg="#0052cc",
        fg="#ffffff",
        command=userinput2_aamod_option_click_enter_button,
        width=8,
    )
    button.grid(row=0, column=1, padx=5, pady=(5, 30))

    window_size_centered(window, 1067, 600)
    window.iconphoto(False, Teriyaki)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.bind("<Return>", userinput2_aamod_option_click_enter_button)
    window.bind("<Escape>", click_close_button)

    window.mainloop()
    return aa_choice, finalchoice


def userinput2_famod(
    aa_options, fa_options, lnk_options, sp_options, aamod_finalchoice, path_to_teriyaki
):
    def click_close_button(event=None):
        nonlocal aa_choice
        aa_choice = ""
        window.destroy()

    def userinput2_famod_click_enter_button(event=None):
        nonlocal aa_choice, fa_choice, lnk_num, lnk_choice, sp_num, sp_choice
        aa_choice = aa_combobox.get()
        fa_choice = fa_combobox.get()
        lnk_num = lnk_spinbox.get()
        lnk_choice = lnk_combobox.get()
        sp_num = sp_spinbox.get()
        sp_choice = sp_combobox.get()

        if aa_choice not in aa_options or aa_choice == "":
            tk.messagebox.showwarning(message="Please choose a valid amino acid entry")
        elif fa_choice not in fa_options or fa_choice == "":
            tk.messagebox.showwarning(message="Please choose a valid fatty acid entry")
        elif lnk_choice not in lnk_options or lnk_choice == "":
            tk.messagebox.showwarning(message="Please choose a valid linker entry")
        elif sp_choice not in sp_options or sp_choice == "":
            tk.messagebox.showwarning(message="Please choose a valid spacer entry")
        else:
            try:
                int(lnk_num)
                int(sp_num)
            except ValueError:
                tk.messagebox.showwarning(
                    message="Please choose a valid linker/spacer number"
                )
            else:
                window.destroy()

    fa_choice = lnk_choice = sp_choice = ""
    lnk_num = sp_num = 0
    aa_choice = ""

    window = tk.Tk()
    window.title("Amino Acid Modification")

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    window_frame = tk.Frame(window)
    window_frame.pack(padx=20, pady=(30, 30))

    aa_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aa_frame.pack(padx=20, pady=(40, 0))

    aa_label1 = tk.Label(aa_frame, text=aamod_finalchoice.upper(), bg="#d9d9d9")
    aa_label1.pack(padx=20, pady=(30, 0))

    aamodaa_entry_frame = tk.LabelFrame(aa_frame, text="", bg="#f0f0f0")
    aamodaa_entry_frame.pack(padx=70, pady=20)

    aachoice_labeltxt = tk.Label(
        aamodaa_entry_frame, text="Choose Amino Acid", bg="#f0f0f0"
    )
    aachoice_labeltxt.pack(padx=10, pady=(15, 0))

    aa_combobox = ttk.Combobox(aamodaa_entry_frame, values=aa_options)
    aa_combobox.set(aa_options[0])
    aa_combobox.pack(padx=20, pady=(10, 20))

    choices_frame = tk.LabelFrame(aa_frame, text="")
    choices_frame.pack(padx=60, pady=(10, 20))

    choices_frame_sub = tk.Frame(choices_frame)
    choices_frame_sub.pack(padx=15)
    choices_frame_sub.columnconfigure(0, minsize=60)
    choices_frame_sub.columnconfigure(1, minsize=40)
    choices_frame_sub.columnconfigure(2, minsize=120)

    fa_label = tk.Label(choices_frame_sub, text="Fatty acid")
    fa_label.grid(row=0, column=0, padx=10, pady=(20, 10), sticky="w")
    fa_combobox = ttk.Combobox(choices_frame_sub, values=fa_options, width=20)
    fa_combobox.set(fa_options[0])
    fa_combobox.grid(row=0, column=1, columnspan=2, padx=5, pady=(20, 10), sticky="w")

    lnk_label = tk.Label(choices_frame_sub, text="Linker")
    lnk_label.grid(row=1, column=0, padx=10, pady=10)
    lnk_spinbox = tk.Spinbox(choices_frame_sub, from_=0, to=5, increment=1, width=4)
    lnk_spinbox.grid(row=1, column=1, padx=(5, 0), pady=10)
    lnk_combobox = ttk.Combobox(choices_frame_sub, values=lnk_options, width=10)
    lnk_combobox.set(lnk_options[0])
    lnk_combobox.grid(row=1, column=2, padx=(0, 5), pady=10)

    sp_label = tk.Label(choices_frame_sub, text="Spacer")
    sp_label.grid(row=2, column=0, padx=10, pady=(0, 20))
    sp_spinbox = tk.Spinbox(choices_frame_sub, from_=0, to=5, increment=1, width=4)
    sp_spinbox.grid(row=2, column=1, padx=(5, 0), pady=(0, 20))
    sp_combobox = ttk.Combobox(choices_frame_sub, values=sp_options, width=10)
    sp_combobox.set(sp_options[0])
    sp_combobox.grid(row=2, column=2, padx=(0, 5), pady=(0, 20))

    button_frame = tk.Frame(aa_frame, bg="#d9d9d9")
    button_frame.pack()

    button_back = tk.Button(
        button_frame, text="Back", bg="#f0f0f0", command=click_close_button, width=8
    )
    button_back.grid(row=0, column=0, padx=5, pady=(15, 30))

    button = tk.Button(
        button_frame,
        text="START",
        bg="#0052cc",
        fg="#ffffff",
        command=userinput2_famod_click_enter_button,
        width=8,
    )
    button.grid(row=0, column=1, padx=5, pady=(15, 30))

    window_size_centered(window, 1067, 600)
    window.iconphoto(False, Teriyaki)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.bind("<Return>", userinput2_famod_click_enter_button)
    window.bind("<Escape>", click_close_button)

    window.mainloop()

    return aa_choice, fa_choice, lnk_num, lnk_choice, sp_num, sp_choice


## USERINPUT 3


def userinput3(peptide_memory, window_status, path_to_teriyaki, path_to_teriyaki_easteregg):

    def click_close_button(event=None):
        nonlocal window_status
        done_input = messagebox.askquestion(
            " ",
            "Are you sure you want to close the program? \n \nATTENTION: Self-added amino acids will be lost! \nSave custom amino acids before closing the program.",
        )
        if done_input == "yes":
            window_status = 4
            main_window.destroy()

    def userinput3_click_startagain_button(event=None):
        nonlocal window_status
        window_status = 1
        main_window.destroy()

    def userinput3_click_back_button(event=None):
        nonlocal back, window_status
        back = True
        window_status = 2
        main_window.destroy()

    def click_getsmiles_button():
        nonlocal molecule
        smiles = Chem.MolToSmiles(molecule)
        pyperclip.copy(smiles)
        messagebox.showinfo("Success", "Value copied to clipboard")

    def click_getsummary_properties():
        nonlocal df_summary_properties
        save_file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=(
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ),
        )
        if save_file_path:
            if save_file_path.endswith(".csv"):
                df_summary_properties.to_csv(save_file_path, index=False)
            elif save_file_path.endswith(".xlsx"):
                df_summary_properties.to_excel(save_file_path, index=False)
            messagebox.showinfo("Save", "File saved successfully")

    def click_savejpg_button():
        nonlocal molecule
        save_jpg_path = filedialog.asksaveasfilename(
            defaultextension=".jpg",
            filetypes=[("JPEG files", "*.jpg"), ("All files", "*.*")],
        )
        if save_jpg_path:
            watermark = "Created with PICKAPEP"
            current_peptide_image_save = Draw.MolToImage(molecule, size=(700, 412))
            current_peptide_draw = ImageDraw.Draw(current_peptide_image_save)
            text_font = ImageFont.truetype("arial.ttf", 15)
            text_width, text_height = current_peptide_draw.textsize(
                watermark, font=text_font
            )
            text_margin = 10
            text_position = (
                current_peptide_image.width - text_width - text_margin,
                current_peptide_image.height - text_height - text_margin,
            )
            current_peptide_draw.text(
                text_position, watermark, fill=(179, 179, 179), font=text_font
            )
            current_peptide_image.save(save_jpg_path, dpi=(300, 300))
            messagebox.showinfo("Save", " Image saved successfully")

    def click_get3Dpymol_button():
        nonlocal peptide_memory
        generate_threedim_pymol(peptide_memory)

    def help_button3_click(event=None):
        subwindow_help3 = tk.Toplevel(main_window, bg="#ffffff")
        subwindow_help3.title("PICKAPEP SUPPORT")

        help3_display = tk.Frame(subwindow_help3, bg="#ffffff")
        help3_display.pack()

        Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)
        subwindow_help3.iconphoto(False, Teriyaki)

        img_Help3 = Image.open(path_to_HS_Userinput3)
        img_Help3_width, img_Help3_height = img_Help3.size
        ratio_size = img_Help3_height / img_Help3_width
        target_width = 800
        img_Help3 = img_Help3.resize((target_width, int(target_width * ratio_size)))
        img_Help3_tk = ImageTk.PhotoImage(img_Help3)
        Help3_display = tk.Label(help3_display, image=img_Help3_tk, bg="#ffffff")
        Help3_display.pack(padx=50, pady=50)

        subwindow_help3.mainloop()

    back = False
    curr_pep = peptide_memory[-1]
    molecule = curr_pep.peptide
    df_summary_properties = create_peptide_summary(peptide_memory)

    main_window = tk.Tk()
    main_window.title("PICKAPEP")

    sequence_easteregg = "".join(
        list(peptide_memory[-1].df_sequence["AminoAcid"][1:-1])
    )
    if sequence_easteregg == "TERIYAKI":
        path_to_teriyaki = path_to_teriyaki_easteregg
        txt_close = "EGG-xit"
    else:
        path_to_teriyaki = path_to_teriyaki
        txt_close = "Close"

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    main_window_frame = tk.Frame(main_window)
    main_window_frame.pack()

    current_structure_label_text = "YOUR FINAL PEPTIDE"
    background_colour = "#d9d9d9"

    win1_current_peptide_frame = tk.LabelFrame(
        main_window_frame, text="", bg=background_colour
    )
    win1_current_peptide_frame.pack(padx=10, pady=(60, 40))

    current_structure_label = tk.Label(
        win1_current_peptide_frame,
        text=current_structure_label_text,
        bg=background_colour,
    )
    current_structure_label.pack(padx=10, pady=(15, 0))

    current_peptide_display = tk.LabelFrame(win1_current_peptide_frame, text="")
    current_peptide_display.pack(padx=30, pady=(10, 15))
    rdCoordGen.AddCoords(molecule)
    current_peptide_image = Draw.MolToImage(molecule, size=(760, 448))
    current_peptide_photo = ImageTk.PhotoImage(current_peptide_image)
    current_peptide_photo_display = tk.Label(
        current_peptide_display, image=current_peptide_photo
    )
    current_peptide_photo_display.pack()

    getinfo_frame_sub = tk.Frame(win1_current_peptide_frame, bg=background_colour)
    getinfo_frame_sub.pack(padx=15, pady=0)

    getsmiles_button = tk.Button(
        getinfo_frame_sub,
        text="Get Smiles",
        bg="#949494",
        fg="#ffffff",
        width=20,
        command=lambda: click_getsmiles_button(),
    )
    getsmiles_button.grid(row=0, column=0, padx=20, pady=(10, 20))
    savejpg_button = tk.Button(
        getinfo_frame_sub,
        text="Save Image",
        bg="#949494",
        fg="#ffffff",
        width=20,
        command=lambda: click_savejpg_button(),
    )
    savejpg_button.grid(row=0, column=1, padx=20, pady=(10, 20))
    getinfo_button = tk.Button(
        getinfo_frame_sub,
        text="Peptide Summary",
        bg="#949494",
        fg="#ffffff",
        width=20,
        command=lambda: click_getsummary_properties(),
    )
    getinfo_button.grid(row=0, column=2, padx=20, pady=(10, 20))
    get3Dpymol_button = tk.Button(
        getinfo_frame_sub,
        text="PyMOL 3D View",
        bg="#949494",
        fg="#ffffff",
        width=20,
        command=lambda: click_get3Dpymol_button(),
    )
    get3Dpymol_button.grid(row=0, column=3, padx=20, pady=(10, 20))

    done_frame = tk.Frame(main_window_frame)
    done_frame.pack()
    startagain_button = tk.Button(
        done_frame,
        text="Start Again",
        bg="#d9d9d9",
        width=12,
        command=lambda: userinput3_click_startagain_button(),
    )
    startagain_button.grid(row=0, column=0)
    back_button = tk.Button(
        done_frame,
        text="Back to Modify",
        bg="#d9d9d9",
        width=12,
        command=lambda: userinput3_click_back_button(),
    )
    back_button.grid(row=0, column=1)
    done_button = tk.Button(
        done_frame,
        text=txt_close,
        bg="#d9d9d9",
        width=12,
        command=lambda: click_close_button(),
    )
    done_button.grid(row=0, column=2)

    icon_frame = tk.Frame(main_window_frame)
    icon_frame.pack(pady=5)
    teriyaki_icon = tk.Label(icon_frame, image=Teriyaki)
    teriyaki_icon.grid(row=0, column=0)
    help_png = Image.open(path_to_helpbutton)
    resize_help_png = help_png.resize((25, 25))
    help_button_png = ImageTk.PhotoImage(resize_help_png)
    help_button = tk.Button(
        icon_frame,
        image=help_button_png,
        command=help_button3_click,
        bd=0,
        highlightthickness=0,
    )
    help_button.grid(row=0, column=1)

    for widget in done_frame.winfo_children():
        widget.grid_configure(padx=5, pady=(0, 15))

    window_size_centered(main_window, 1422, 800)
    main_window.iconphoto(False, Teriyaki)

    main_window.protocol("WM_DELETE_WINDOW", click_close_button)

    main_window.bind("<Return>", click_close_button)
    main_window.bind("<Escape>", click_close_button)
    main_window.bind("<BackSpace>", userinput3_click_back_button)
    main_window.bind("<Delete>", userinput3_click_startagain_button)
    main_window.bind("<question>", help_button3_click)

    main_window.mainloop()

    return peptide_memory, back, window_status


## ADDITIONAL ##


def input_singlechoice(choices):

    def click_close_button(event=None):
        nonlocal finalchoice
        finalchoice = ""
        window.destroy()

    def save_choice(event=None):
        nonlocal finalchoice
        finalchoice = selected_choice.get()
        window.destroy()

    aa_choice = ""

    window = tk.Tk()
    window.title("Cyclization Information")

    Teriyaki = ImageTk.PhotoImage(file=path_to_teriyaki)

    window_frame = tk.Frame(window)
    window_frame.pack(padx=20, pady=(30, 30))

    aa_frame = tk.LabelFrame(window_frame, text="", bg="#d9d9d9")
    aa_frame.pack(padx=20, pady=(60, 0))

    aa_label1 = tk.Label(
        aa_frame, text="Which reaction should be performed?", bg="#d9d9d9"
    )
    aa_label1.pack(padx=20, pady=(30, 0))

    aamodaa_entry_frame = tk.LabelFrame(aa_frame, text="", bg="#f0f0f0")
    aamodaa_entry_frame.pack(padx=70, pady=15)

    selected_choice = tk.StringVar(value=choices[0])
    finalchoice = []

    for idx, choice in enumerate(choices):
        if idx == 0:
            y_value = (20, 5)
        elif idx == len(choices) - 1:
            y_value = (5, 20)
        else:
            y_value = 5
        radiobutton = tk.Radiobutton(
            aamodaa_entry_frame, text=choice, variable=selected_choice, value=choice
        )
        radiobutton.grid(row=idx, column=0, sticky="W", padx=20, pady=y_value)

    button_frame = tk.Frame(aa_frame, bg="#d9d9d9")
    button_frame.pack()

    button_back = tk.Button(
        button_frame, text="Back", bg="#f0f0f0", command=click_close_button, width=8
    )
    button_back.grid(row=0, column=0, padx=5, pady=(15, 30))

    button = tk.Button(
        button_frame,
        text="START",
        bg="#0052cc",
        fg="#ffffff",
        command=save_choice,
        width=8,
    )
    button.grid(row=0, column=1, padx=5, pady=(15, 30))

    window_size_centered(window, 1067, 600)
    window.iconphoto(False, Teriyaki)

    window.protocol("WM_DELETE_WINDOW", click_close_button)

    window.bind("<Return>", save_choice)
    window.bind("<Escape>", click_close_button)

    window.mainloop()

    return finalchoice


# Userinput 3 when fast calculation is done
def userinput3_fastcalc(calc_table, calc_table_path, skip):
    count_done_peptides = (calc_table["Status"] == "Done").sum()
    count_empty_peptides = (calc_table["Status"] == "").sum()
    done_peptides = count_done_peptides - skip
    tocalc_peptides = len(calc_table) - skip - count_empty_peptides
    str_message = (
        "Calculations done: "
        + str(done_peptides)
        + "/"
        + str(tocalc_peptides)
        + " Peptides calculated successfully. \nDo you want to save the file?"
    )
    save_result = messagebox.askquestion("", str_message)
    add_text = ""
    if save_result == "yes":
        calc_table_file_name = os.path.splitext(os.path.basename(calc_table_path))[0]
        new_calc_table_file_name = f"{calc_table_file_name}_PICKAPEP"
        new_calc_table_file_path = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(calc_table_path),
            initialfile=new_calc_table_file_name,
            defaultextension=".xlsx",
            filetypes=(
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ),
        )
        if new_calc_table_file_path:
            if new_calc_table_file_path.endswith(".csv"):
                calc_table.to_csv(new_calc_table_file_path, index=False)
            elif new_calc_table_file_path.endswith(".xlsx"):
                calc_table.to_excel(new_calc_table_file_path, index=False)
            add_text = "File Saved Successfully! "
    again_result = messagebox.askquestion(
        "",
        add_text
        + "Do you want to go back to the start? \n\nATTENTION: Self-added amino acids will be lost! \nSave custom amino acids before closing the program.",
    )
    if again_result == "yes":
        window_status = 1
    else:
        window_status = 4

    return window_status


def keepmapnum(reactant1, reactant2, product, reaction):
    """Reassigns mapnums assigned to reactants before reaction to product.
    Mapnums are lost during reaction in atoms part of the reaction template.

    Parameters
    ----------
    reactant1 : rdkit.Chem.rdchem.Mol
        RDkit mol object which was reactant 1 in reaction
    reactant2 : rdkit.Chem.rdchem.Mol
        RDkit mol object which was reactant 2 in reaction
    product : rdkit.Chem.rdchem.Mol
        RDkit mol object which was product of the reaction
    reaction : rdkit.Chem.rdChemReactions.ChemicalReaction
        RDkit chemical reaction which was performed with reactant1 and reactant2, forming product

    Returns
    ------
    product : rdkit.Chem.rdchem.Mol
        RDkit mol object with reassigned mapnums
    """

    # 1) figure out which reactant each atom came from
    atomMapToReactantMap = {}
    reactants = (reactant1, reactant2)
    for ri in range(reaction.GetNumReactantTemplates()):
        rt = reaction.GetReactantTemplate(ri)
        for atom in rt.GetAtoms():
            if atom.GetAtomMapNum():
                atomMapToReactantMap[atom.GetAtomMapNum()] = ri

    # 2) from that reactant, copy over the mapnums to product
    for atom in product.GetAtoms():
        if atom.HasProp("old_mapno"):
            reactant_this_came_from = reactants[
                atomMapToReactantMap[int(atom.GetProp("old_mapno"))]
            ]
            fitting_atoms = [
                ra
                for ra in reactant_this_came_from.GetAtoms()
                if ra.GetIdx() == atom.GetIntProp("react_atom_idx")
            ]  # todo we could use a dict
            assert len(fitting_atoms) == 1
            orig_atom = fitting_atoms[0]
            assert atom.GetSymbol() == orig_atom.GetSymbol()  # sanity check
            atom.SetIntProp(
                "mapnum", orig_atom.GetIntProp("mapnum")
            )  # copy over from old atom

    return product


def find_mapnum_product(output_tuple, mapnumprop):
    """Find in output_tuple products, with the correct reaction sites.

    Parameters
    ----------
    output_tuple : tuple
        tuple containing reaction products resulting from an RDkit chemical reaction
    mapnumprop : list
        list of mapnums of atoms at desired reaction site

    Returns
    ------
    products : tuple
        tuple (same format as output_tuple) containing reaction products with correct reaction sites
    """

    rx_output = tuple(element[0] for element in output_tuple)
    product_store = []
    for i in rx_output:
        mapnum_atoms_idx = [
            atom.GetIdx()
            for atom in i.GetAtoms()
            if atom.HasProp("mapnum") and atom.GetIntProp("mapnum") in mapnumprop
        ]
        rx_atoms_idx = [
            atom.GetIdx() for atom in i.GetAtoms() if atom.HasProp("old_mapno")
        ]
        if all(value in mapnum_atoms_idx for value in rx_atoms_idx):
            product_store = product_store + [i]
    products = tuple((value,) for value in product_store)

    return products


def num_equiv_solutions(output_tuple):
    """Determine number of equivalent products by an RDkit chemical reaction.

    Parameters
    ----------
    output_tuple : tuple
        tuple containing reaction products resulting from an RDkit chemical reaction

    Returns
    ------
    num_solutions : int
        number of equivalent solutions generated by RDkit chemical reaction
    """

    rx_output = tuple(element[0] for element in output_tuple)
    molcleaned_products = [(Chem.MolToSmiles(x)) for x in rx_output]
    num_solutions = len(set(molcleaned_products))

    return num_solutions


def remove_intraaacycl_product(output_tuple, aa_mapnumprop):
    """Cleans up cyclization products: removes products where not both reactions sites chosen for reaction were involved in reaction. (Can happen when amino acid has two reaction sites e.g. B5).

    Parameters
    ----------
    output_tuple : tuple
        tuple containing reaction products resulting from an RDkit chemical reaction
    aa_mapnumprop : list
        list of mapnums of each reaction site atoms at desired reaction sites

    Returns
    ------
    products : tuple
        tuple (same format as output_tuple) containing reaction products where all chosen reaction sites were involved in reaction
    """

    rx_output = tuple(element[0] for element in output_tuple)
    product_store = []
    for i in rx_output:
        mapnum_atoms_idx = [
            [
                atom.GetIdx()
                for atom in i.GetAtoms()
                if atom.HasProp("mapnum") and atom.GetIntProp("mapnum") in sublist
            ]
            for sublist in aa_mapnumprop
        ]
        rx_atoms_idx = [
            atom.GetIdx() for atom in i.GetAtoms() if atom.HasProp("old_mapno")
        ]
        if all(
            any(element in rx_atoms_idx for element in sublist)
            for sublist in mapnum_atoms_idx
        ):
            product_store = product_store + [i]
    products = tuple((value,) for value in product_store)
    return products


def assign_pepinfo(peptide_mol, df_sequence):
    """Calculates peptide parameters and saves them in a DataFrame

    Parameters
    ----------
    peptide_mol : rdkit.Chem.rdchem.Mol
        RDkit mol object of the peptide
    df_sequence : pd.DataFrame
        DataFrame with peptide sequence information

    Returns
    ------
    df_peptide_info: pd.DataFrame
        DataFrame containing calculated peptide parameters
    """

    chemical_formular = AllChem.CalcMolFormula(peptide_mol)
    molwt = Chem.Descriptors.MolWt(peptide_mol)
    exactmolwt = Chem.Descriptors.ExactMolWt(peptide_mol)
    exactsinglecharge = exactmolwt + Chem.Descriptors.ExactMolWt(
        Chem.MolFromSmiles("[H]")
    )
    exactdoublecharge = (
        exactmolwt + 2 * Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles("[H]"))
    ) / 2
    exacttriplecharge = (
        exactmolwt + 3 * Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles("[H]"))
    ) / 3
    exactquadruplecharge = (
        exactmolwt + 4 * Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles("[H]"))
    ) / 4
    logp = Chem.Descriptors.MolLogP(peptide_mol)
    TPSA_val = Chem.Descriptors.TPSA(peptide_mol)
    rotbonds_num = Chem.Descriptors.NumRotatableBonds(peptide_mol)
    HBA_num = Chem.Descriptors.NumHAcceptors(peptide_mol)
    HBD_num = Chem.Descriptors.NumHDonors(peptide_mol)
    Arom_cont = Chem.Descriptors.NumAromaticRings(peptide_mol)
    aa_num = len(df_sequence) - 2  # remove N- and C-terminus

    peptide_info = {
        "Info": [
            "Chemical Formula",
            "Molecular Weight",
            "Exact Weight",
            "[M+H]+",
            "[M+2H]2+",
            "[M+3H]3+",
            "[M+4H]4+",
            "Number of Amino Acids",
            "LogP",
            "TPSA",
            "Rotatable Bonds",
            "Number of HBA",
            "Number of HBD",
            "Aromatic Rings",
        ],
        "Value": [
            chemical_formular,
            str(round(molwt, 2)),
            str(round(exactmolwt, 4)),
            str(round(exactsinglecharge, 4)),
            str(round(exactdoublecharge, 4)),
            str(round(exacttriplecharge, 4)),
            str(round(exactquadruplecharge, 4)),
            str(aa_num),
            str(round(logp, 2)),
            str(round(TPSA_val, 2)),
            str(int(rotbonds_num)),
            str(int(HBA_num)),
            str(int(HBD_num)),
            str(int(Arom_cont)),
        ],
    }

    df_peptide_info = pd.DataFrame(peptide_info)

    return df_peptide_info


def check_sequence_input(sequence, allowed, df_aminoacids):
    """Checks entered sequence, removes not allowed symbols, sorts input into amino acids (according to capital letters).
    Checks if amino acids are available in library and assigns smiles code. If amino acid is not available in library, an errormessage containing the unknown amino acid is returned.

    Parameters
    ----------
    sequence : str
        Userinput amino acid sequence
    allowed : str
        Regular expression of symbols allowed as elements in the sequence (str), contains upper- and lowercase letters and numbers
    df_aminoacids: pd.DataFrame
        DataFrame of amino acid library linking entry code with smiles code

    Returns
    ------
    aaseq : list
        list of amino acid entry codes
    aa_smiles : list
        list of smiles code of amino acids to assemble
    errormessage : str
        'none' if all amino acids are found in library, otherwise it returns the string of the unknown amino acid to inform the user
    """

    errormessage = "none"

    sequence = re.sub(allowed, "", str(sequence))
    seq_ascii = [ord(i) for i in sequence]

    join_elem = []
    seq = []

    # sort input according to capital letters (-> marks new amino acid ASCII 65-90 A-Z)
    # assign small letters and numbers to capital letter (ASCII 97-122 a-z, ASCII 48-57 0-9)

    for i in range(len(seq_ascii)):
        elem = seq_ascii[i]
        if elem >= 65 and elem <= 90:
            join_elem = [elem]
        elif (elem >= 97 and elem <= 122) or (elem >= 48 and elem <= 57):
            join_elem.append(elem)
        if i + 1 < len(seq_ascii):
            if len(join_elem) != 0:
                if seq_ascii[i + 1] >= 65 and seq_ascii[i + 1] <= 90:
                    seq.append(join_elem)
        else:
            seq.append(join_elem)

    # Convert amino acid sequence input sorted back to text
    aaseq = []
    for aa in seq:
        aatxt = "".join(chr(i) for i in aa)
        aaseq.append(aatxt)

    # Check sequence for available aminoacids, unify all to oneletter notation and assign smiles code
    aa_smiles = [0] * len(aaseq)

    for i in range(len(aaseq)):
        aa = aaseq[i]
        if any(string == aa for string in df_aminoacids.loc[:, "Entry_Letter"]):
            idx = pd.Index(df_aminoacids.loc[:, "Entry_Letter"]).get_loc(aa)
            aa_smiles[i] = df_aminoacids.loc[idx, "Smiles"]
        elif any(string == aa for string in df_aminoacids.loc[:, "Entry_Code"]):
            idx = pd.Index(df_aminoacids.loc[:, "Entry_Code"]).get_loc(aa)
            aaseq[i] = df_aminoacids.loc[idx, "Entry_Letter"]
            aa_smiles[i] = df_aminoacids.loc[idx, "Smiles"]
        else:
            errormessage = f"Unknown aminoacid: {aa}"

    return aaseq, aa_smiles, errormessage


def create_linear_peptide(aaseq, aa_smiles, df_sequence_overview, errormessage):
    """Assembles amino acid building blocks to linear peptide and stores assigned mapnum for each amino acid building block to ensure regioselectivity for further modifications.
    Checks generated molecule after each amino acid attachement and returns errormessage if linear peptide assembly did not work.

    Parameters
    ----------
    aaseq : list
        list of amino acid entry codes
    aa_smiles : list
        list of smiles code of amino acids
    df_sequence_overview: pd.DataFrame
        DataFrame to store sequence information including mapnum
    errormessage : str
        current errormessage status to allow update if needed

    Returns
    ------
    chain : rdkit.Chem.rdchem.Mol
        RDkit mol object of the linear peptide
    df_sequence_overview : pd.DataFrame
        DataFrame updated with linear sequence information
    mapnum_store : list
        list of mapnums stored along peptide assembly
    highest_mapnum : int
        current highest mapnum assigned to an atom
    errormessage :
        'none' if all amino acids are successfully attached, otherwise it returns 'Linear peptide assembly failed' to inform the user
    """

    aa_mols = [Chem.MolFromSmiles(i) for i in aa_smiles]

    # Peptide bond formation
    smarts_reaction_peptidebond = "[NX3:1][CX4:2][CX3:3](=[O:4])[O].[NX3:5][CX4:6][CX3:7](=[O:8])[O:9]>>[N:1][C:2][C:3](=[O:4])[N:5][C:6][C:7](=[O:8])[O:9]"
    reaction_peptidebond = rdChemReactions.ReactionFromSmarts(
        smarts_reaction_peptidebond
    )

    chain = aa_mols[0]

    # Create mapnum so amino acids can be distiguished
    [a.SetIntProp("mapnum", i + 1) for i, a in enumerate(chain.GetAtoms())]
    highest_mapnum = chain.GetNumAtoms()
    mapnum_store = [0, highest_mapnum]

    for i in range(1, len(aa_mols)):
        # assign mapnums to all atoms in the next molecule we're adding
        aa_mol = aa_mols[i]
        [a.ClearProp("old_mapno") for a in chain.GetAtoms()]
        [
            a.SetIntProp("mapnum", highest_mapnum + i + 1)
            for i, a in enumerate(aa_mol.GetAtoms())
        ]

        # Amino acid attachement reaction, clean output, store mapnum
        chain_rx = reaction_peptidebond.RunReactants((chain, aa_mol))
        chain_rx = tuple(
            (value,)
            for value in [
                keepmapnum(chain, aa_mol, solution[0], reaction_peptidebond)
                for solution in chain_rx
            ]
        )

        try:
            if num_equiv_solutions(chain_rx) == 1:
                chain_ext = chain_rx[0][0]
                Chem.SanitizeMol(chain_ext)
                chain = chain_ext
            else:
                errormessage = (
                    "Linear peptide assembly failed at amino acid "
                    + str(i)
                    + " "
                    + aaseq[i]
                )
                return (
                    chain,
                    df_sequence_overview,
                    mapnum_store,
                    highest_mapnum,
                    errormessage,
                )
        except:
            errormessage = (
                "Linear peptide assembly failed at amino acid "
                + str(i)
                + " "
                + aaseq[i]
            )
            return (
                chain,
                df_sequence_overview,
                mapnum_store,
                highest_mapnum,
                errormessage,
            )

        # Update mapnum if reaction worked
        highest_mapnum += aa_mol.GetNumAtoms()
        mapnum_store.append(highest_mapnum)

    # update df_sequence_overview with linear peptide information
    df_sequence_overview.loc[0, "AminoAcid"] = "N-Term"

    for i in range(1, len(aaseq) + 1):
        aa_mapnum = list(range(mapnum_store[i - 1] + 1, mapnum_store[i] + 1))
        df_sequence_overview.loc[i, ["AminoAcid", "Mapnum"]] = aaseq[i - 1], aa_mapnum

    df_sequence_overview.loc[len(aaseq) + 1, "AminoAcid"] = "C-Term"

    for index, value in enumerate(df_sequence_overview["AminoAcid"]):
        if index == 0 or index == len(df_sequence_overview) - 1:
            df_sequence_overview.at[index, "AminoAcidInfo"] = f" {value} "
        else:
            df_sequence_overview.at[index, "AminoAcidInfo"] = f" {index} {value} "

    return chain, df_sequence_overview, mapnum_store, highest_mapnum, errormessage


def attach_nterminus(
    chain,
    df_sequence_overview,
    mapnum_store,
    highest_mapnum,
    dict_nterm,
    nterm,
    errormessage,
):
    """N-terminus modification of linear peptide.

    Parameters
    ----------
    chain : rdkit.Chem.rdchem.Mol
        RDkit mol object of the linear peptide
    df_sequence_overview: pd.DataFrame
        DataFrame storing linear sequence information including mapnum
    mapnum_store : list
        list of mapnums stored along peptide assembly
    highest_mapnum : int
        current highest mapnum assigned to an atom
    dict_nterm : dict
        dictionary containing available n termini
    nterm : str
        userinput nterm
    errormessage : str
        current errormessage status to allow update if needed

    Returns
    ------
    chain : rdkit.Chem.rdchem.Mol
        RDkit mol object of the linear peptide with nterminus
    df_sequence_overview : pd.DataFrame
        DataFrame updated with N-terminus information
    mapnum_store : list
        list of mapnums stored including N-terminus attached
    highest_mapnum : int
        current highest mapnum assigned to an atom
    errormessage :
        'none' if N-terminus attachement succesful, otherwise it returns 'N-terminus attachement failed' to inform the user
    """

    tail_mol = Chem.MolFromSmiles(dict_nterm[nterm]["Smiles"])
    rx_nterm = rdChemReactions.ReactionFromSmarts(dict_nterm[nterm]["Reaction_Smarts"])
    [a.ClearProp("old_mapno") for a in chain.GetAtoms()]
    [
        a.SetIntProp("mapnum", highest_mapnum + j + 1)
        for j, a in enumerate(tail_mol.GetAtoms())
    ]

    # Store right mymap num in tail
    highest_mapnum += tail_mol.GetNumAtoms()
    mapnum_store.append(highest_mapnum)
    aa_mapnum = list(range(mapnum_store[-2] + 1, mapnum_store[-1] + 1))
    df_sequence_overview.at[0, "Mapnum"] = aa_mapnum
    df_sequence_overview.loc[0, "AminoAcid"] = nterm
    df_sequence_overview.at[0, "AminoAcidInfo"] = (
        df_sequence_overview.loc[0, "AminoAcidInfo"] + " " + nterm
    )

    # Make reaction for N term Modification
    chain_nterm = rx_nterm.RunReactants((chain, tail_mol))
    chain_nterm = tuple(
        (value,)
        for value in [
            keepmapnum(chain, tail_mol, solution[0], rx_nterm)
            for solution in chain_nterm
        ]
    )

    # Check for right product
    rx_nterm_mapnum = (
        df_sequence_overview.loc[0, "Mapnum"] + df_sequence_overview.loc[1, "Mapnum"]
    )
    chain_nterm = find_mapnum_product(chain_nterm, rx_nterm_mapnum)

    if num_equiv_solutions(chain_nterm) == 1:
        chain_ext = chain_nterm[0][0]
        Chem.SanitizeMol(chain_ext)
        chain = chain_ext
    else:
        chain = ""
        errormessage = "N-terminus attachement failed"
    return chain, df_sequence_overview, mapnum_store, highest_mapnum, errormessage


def attach_cterminus(
    chain,
    df_sequence_overview,
    mapnum_store,
    highest_mapnum,
    dict_cterm,
    cterm,
    errormessage,
):
    """N-terminus modification of linear peptide.

    Parameters
    ----------
    chain : rdkit.Chem.rdchem.Mol
        RDkit mol object of the linear peptide
    df_sequence_overview: pd.DataFrame
        DataFrame storing linear sequence information including mapnum
    mapnum_store : list
        list of mapnums stored along peptide assembly
    highest_mapnum : int
        current highest mapnum assigned to an atom
    dict_cterm : dict
        dictionary containing available c termini
    cterm : str
        userinput cterm
    errormessage : str
        current errormessage status to allow update if needed

    Returns
    ------
    chain : rdkit.Chem.rdchem.Mol
        RDkit mol object of the linear peptide with nterminus
    df_sequence_overview : pd.DataFrame
        DataFrame updated with C-terminus information
    mapnum_store : list
        list of mapnums stored including C-terminus attached
    highest_mapnum : int
        current highest mapnum assigned to an atom
    errormessage :
        'none' if C-terminus attachement succesful, otherwise it returns 'C-terminus attachement failed' to inform the user
    """

    cterm_count = len(df_sequence_overview) - 1

    head_mol = Chem.MolFromSmiles(dict_cterm[cterm]["Smiles"])
    rx_cterm = rdChemReactions.ReactionFromSmarts(dict_cterm[cterm]["Reaction_Smarts"])
    [a.ClearProp("old_mapno") for a in chain.GetAtoms()]
    [
        a.SetIntProp("mapnum", highest_mapnum + j + 1)
        for j, a in enumerate(head_mol.GetAtoms())
    ]

    # Store right mymap num in head
    highest_mapnum += head_mol.GetNumAtoms()
    mapnum_store.append(highest_mapnum)
    aa_mapnum = list(range(mapnum_store[-2] + 1, mapnum_store[-1] + 1))
    df_sequence_overview.at[cterm_count, "Mapnum"] = aa_mapnum
    df_sequence_overview.loc[cterm_count, "AminoAcid"] = cterm
    df_sequence_overview.at[cterm_count, "AminoAcidInfo"] = (
        df_sequence_overview.loc[cterm_count, "AminoAcidInfo"] + " " + cterm
    )

    # Make reaction for C term Modification
    chain_cterm = rx_cterm.RunReactants((chain, head_mol))
    chain_cterm = tuple(
        (value,)
        for value in [
            keepmapnum(chain, head_mol, solution[0], rx_cterm)
            for solution in chain_cterm
        ]
    )

    # Check for right product
    rx_cterm_mapnum = (
        df_sequence_overview.loc[cterm_count - 1, "Mapnum"]
        + df_sequence_overview.loc[cterm_count, "Mapnum"]
    )
    chain_cterm = find_mapnum_product(chain_cterm, rx_cterm_mapnum)

    if num_equiv_solutions(chain_cterm):
        chain_ext = chain_cterm[0][0]
        Chem.SanitizeMol(chain_ext)
        chain = chain_ext
    else:
        chain = ""
        errormessage = "C-terminus attachement failed"
    return chain, df_sequence_overview, mapnum_store, highest_mapnum, errormessage


def sort_cycl_htinput(cycl_htinput, df_sequence_overview):
    """Sorts the string input for cyclization for a given peptide in the high throughput calculation mode putting it into a form PICKAPEP can use further.
    Avoids wrong input processing.

    Parameters
    ----------
    cycl_htinput : str
        Cyclization input for a peptide in high throughput calculation mode.
    df_sequence_overview : pd.DataFrame
        DataFrame storing linear sequence information including mapnum

    Returns
    ------
    cycl_sorted_list : list
        List storing the necessary information for peptide cyclization (position numbers of involved amino acids and specification) for as many cyclization reactions as defined.
    """

    cycl_sorted_list = []
    cycl_htinput = re.sub(r"[^a-zA-Z0-9,;]", "", str(cycl_htinput))
    cycl_htinput = cycl_htinput.lower()

    if cycl_htinput and not cycl_htinput.lower() == "none":
        cycl_rx_info = cycl_htinput.split(";")
        cycl_rx_info = [item for item in cycl_rx_info if item != ""]

        for elem in cycl_rx_info:
            cycl_rx_var = elem.split(",")
            cycl_rx_var = [item for item in cycl_rx_var if item != ""]
            numbers = []
            name = ""
            for var in cycl_rx_var:
                if var.lower() == "head":
                    num = len(df_sequence_overview) - 1
                    numbers.append(num)
                elif var.lower() == "tail":
                    num = 0
                    numbers.append(num)
                elif var.isdigit():
                    num = int(var)
                    if num <= len(df_sequence_overview):
                        numbers.append(num)
                else:
                    name = var
            cycl_sorted_list.append([numbers, name])
    return cycl_sorted_list


def sort_aamod_htinput(aamod_htinput, df_sequence_overview, df_modifications):
    """Sorts the string input for modifications for a given peptide in the high throughput calculation mode putting it into a form PICKAPEP can use further.
    Avoids wrong input processing.

    Parameters
    ----------
    aamod_htinput : str
        Modification input for a peptide in high throughput calculation mode.
    df_sequence_overview : pd.DataFrame
        DataFrame storing linear sequence information including mapnum.
    df_modifications : pandas.core.frame.Dataframe
        Contains all possible modification reactions.


    Returns
    ------
    aamod_sorted_list : list
        List storing the necessary information for peptide modifications (position numbers of involved amino acids and specification) for as many modification reactions as defined.
    """

    aamod_sorted_list = []
    aamod_htinput = re.sub(r"[^a-zA-Z0-9,;]", "", str(aamod_htinput))
    aamod_htinput = aamod_htinput.lower()

    if aamod_htinput and not aamod_htinput.lower() == "none":
        aamod_rx_info = aamod_htinput.split(";")
        aamod_rx_info = [item for item in aamod_rx_info if item != ""]

        for elem in aamod_rx_info:
            aamod_rx_var = elem.split(",")
            aamod_rx_var = [item for item in aamod_rx_var if item != ""]
            numbers = []
            fa_ht_info = []
            name = ""

            if (
                aamod_rx_var[1] == "fa"
                or aamod_rx_var[0] == "fa"
                or any(
                    (
                        not df_modifications[
                            df_modifications["Variable"] == str(fa_var)
                        ].empty
                    )
                    and (
                        df_modifications[df_modifications["Variable"] == str(fa_var)][
                            "Type"
                        ].iloc[0]
                        == "fattyac"
                    )
                    for fa_var in aamod_rx_var[:2]
                )
            ) and len(aamod_rx_var) > 2:
                for var in aamod_rx_var[:2]:
                    if var.isdigit():
                        num = int(var)
                        numbers = [num]
                    else:
                        if name == "":
                            name = var

                for index in range(2, len(aamod_rx_var)):
                    var = aamod_rx_var[index]
                    if var.isdigit():
                        continue
                    else:
                        fa_sub_name = var
                        if index != 2 and aamod_rx_var[index - 1].isdigit():
                            fa_sub_num = int(aamod_rx_var[index - 1])
                        else:
                            fa_sub_num = 1
                    fa_ht_info.append([fa_sub_num, fa_sub_name])
                aamod_sorted_list.append([numbers, name, fa_ht_info])
            else:
                for var in aamod_rx_var:
                    if var.isdigit():
                        num = int(var)
                        numbers = [num]
                    else:
                        if name == "":
                            name = var
                aamod_sorted_list.append([numbers, name])

    return aamod_sorted_list


def sort_cycl_input(
    df_sequence_overview,
    aa_count,
    cyclization_finalchoice,
    df_reactions,
    userinput2_cycl,
):
    """Sorts userinput for cyclization reaction and calls for further userinput if needed.

    Parameters
    ----------
    df_sequence_overview : pd.DataFrame
        DataFrame storing linear sequence information including mapnum.
    aa_count : int
        Number of amino acids in peptide.
    cyclization_finalchoice: str
        Chosen cyclization overcategory (Head to Tail, Head/Tail to Sidechain, Sidechain to Sidechain, Two component, Bicycle)
    df_reactions: pandas.core.frame.DataFrame
        DataFrame storing all available cyclization reactions.
    userinput2_cycl: function
        Asks for further userinput depending on the chosen cyclization type if needed.


    Returns
    ------
    aa_choice_output: list
        List of strings of chosen amino acids for cyclization.
    aa_choice_num_output: list
        List of integers containing the aminoacid numbers of chosen amino acids for cyclization.
    additional_needed: bool
        Boolean containing the information if additional userinput to define the cyclization reaction besides the amino acid numbers is needed.
    additional_info: str
        String of the additional info needed to ask for in the userinput.
    additional_output: str
        String containing the additional info needed for the cyclization reaction given by userinput.
    skip_cycl: bool
        Boolean containing the information if the reaction should be performed or not (in case the reaction is not possible or user decided to not give needed input). If false, the user gets back to the previous peptide.
    """

    # Variables that should be reset but could be handy for cycl as well as modification
    additional_needed = False
    input_needed = True
    additional_options = []
    additional_info = ""
    additional_output = ""

    aa_options = df_sequence_overview.loc[1:aa_count, "AminoAcidInfo"].tolist()
    aa_head = [df_sequence_overview.loc[aa_count + 1, "AminoAcidInfo"]]
    aa_tail = [df_sequence_overview.loc[0, "AminoAcidInfo"]]

    # Set all variables for userinput depending on cyclization type
    if cyclization_finalchoice == "Head to Tail":
        input_needed = False
        aa_choice_output = [aa_tail[0], aa_head[0]]
        aa_choice_num_output = [0, aa_count + 1]

    elif cyclization_finalchoice == "Head to Sidechain":
        aa_options_input = [aa_options, aa_head]

    elif cyclization_finalchoice == "Tail to Sidechain":
        aa_options_input = [aa_tail, aa_options]

    elif cyclization_finalchoice == "Sidechain to Sidechain":
        aa_options_input = [aa_options, aa_options]

    elif cyclization_finalchoice == "Two Component":
        aa_options_input = [aa_options, aa_options]
        additional_needed = True
        additional_options = df_reactions.loc[
            df_reactions["Type"] == "twocomp", "Name"
        ].tolist()
        additional_info = "Second Component for Two Component Stapling"

    elif cyclization_finalchoice == "Bicycle Formation":
        aa_options_input = [aa_options, aa_options, aa_options]
        additional_needed = True
        additional_options = df_reactions.loc[
            df_reactions["Type"] == "bicycle", "Name"
        ].tolist()
        additional_info = "Bicycle Formation"

    skip_cycl = False
    if input_needed:
        aa_choice_output, additional_output, skip_cycl = userinput2_cycl(
            cyclization_finalchoice,
            aa_options_input,
            additional_needed,
            additional_options,
            additional_info,
            path_to_teriyaki,
        )
        aa_choice_num_output = []
        if skip_cycl:
            return (
                aa_choice_output,
                aa_choice_num_output,
                additional_needed,
                additional_info,
                additional_output,
                skip_cycl,
            )
        for elem in aa_choice_output:
            if elem[0:7] == " N-Term":
                choice_num = 0
            elif elem[0:7] == " C-Term":
                choice_num = aa_count + 1
            else:
                choice_num = int([part for part in elem.split(" ") if part][0])
            aa_choice_num_output.append(choice_num)

    return (
        aa_choice_output,
        aa_choice_num_output,
        additional_needed,
        additional_info,
        additional_output,
        skip_cycl,
    )


def get_possible_cycl_rx(
    aa_choice_num_output,
    df_sequence_overview,
    df_reactions,
    chain,
    cyclization_finalchoice,
    rx_finalchoice,
    cycl_choice_var,
    additional_output,
    calc_type,
):
    """Find cyclization reactions possible with the current peptide between the chosen amino acids (with additional component if given).
    In case of calc_type == fast it checks if the desired reaction is possible between the chosen amino acids.

    Parameters
    ----------
    aa_choice_num_output: list
        List of integers containing the aminoacid numbers of chosen amino acids for cyclization.
    df_sequence_overview : pd.DataFrame
        DataFrame storing linear sequence information including mapnum.
    df_reactions: pandas.core.frame.DataFrame
        DataFrame storing all available cyclization reactions.
    chain: rdkit.Chem.rdchem.Mol
        RDkit mol object of the current peptide
    cyclization_finalchoice: str
        Chosen cyclization overcategory (Head to Tail, Head/Tail to Sidechain, Sidechain to Sidechain, Two component, Bicycle)
    rx_finalchoice: str
        Contain in case of calc_type == fast the final reaction the user intends to perform. In case of calc_type == slow not defined yet.
    cycl_choice_var: str
        variable notation of cyclization finalchoice
    additional_output: str
        String containing the additional info needed for the cyclization reaction given by previous userinput.
    calc_type: str
        Either slow or fast depending on the calculation mode the user has chosen.


    Returns
    ------
    possible_rx: list
        List of reactions possible between the amino acids chosen for the given reaction type.
    rx_mapnum: list
        List of mapnums of the involved amino acids.
    aa_mapnum: list
        List of list with mapnums per amino acid involved in intended cyclization reaction.
    additional_output: str
        String containing the additional info needed given by userinput.
    """

    rx_mapnum = []
    aa_mapnum = []
    for i in aa_choice_num_output:
        rx_mapnum = rx_mapnum + df_sequence_overview.loc[i, "Mapnum"]
        aa_mapnum = aa_mapnum + [df_sequence_overview.loc[i, "Mapnum"]]
    rx_mapnum_atomidx = [
        atom.GetIdx()
        for atom in chain.GetAtoms()
        if atom.GetIntProp("mapnum") in rx_mapnum
    ]

    possible_rx = []
    if (
        cyclization_finalchoice == "Two Component"
        or cyclization_finalchoice == "Bicycle Formation"
    ):
        if calc_type == "fast":
            additional_output = (
                df_reactions[df_reactions["Variable"] == rx_finalchoice]
                .reset_index(drop=True)
                .loc[0, "Name"]
            )
        subdf_reaction = df_reactions[
            df_reactions["Name"] == additional_output
        ].reset_index(drop=True)
        rx_finalchoice = subdf_reaction.iloc[0]["Variable"]
        x = rdChemReactions.ReactionFromSmarts(
            subdf_reaction.iloc[0]["Reaction_Smarts"]
        )
        for match in chain.GetSubstructMatches(x.GetReactantTemplate(0)):
            if all(idx in rx_mapnum_atomidx for idx in match):
                if rx_finalchoice not in possible_rx:
                    possible_rx = possible_rx + [rx_finalchoice]
    else:  # Head to Tail, Side Chain to Side Chain, Head to Side Chain, Tail to Side Chain
        for index, row in df_reactions[
            df_reactions["Type"] == cycl_choice_var
        ].iterrows():
            x = rdChemReactions.ReactionFromSmarts(
                df_reactions[df_reactions["Type"] == cycl_choice_var].loc[index][
                    "Reaction_Smarts"
                ]
            )
            for match in chain.GetSubstructMatches(x.GetReactantTemplate(0)):
                if all(idx in rx_mapnum_atomidx for idx in match):
                    if df_reactions.loc[index, "Variable"] not in possible_rx:
                        possible_rx = possible_rx + [
                            df_reactions.loc[index, "Variable"]
                        ]
    return possible_rx, rx_mapnum, aa_mapnum, additional_output


def perform_cycl_rx(
    df_reactions,
    cycl_choice_var,
    rx_finalchoice,
    chain_reactant,
    additional_reactant,
    highest_mapnum,
    mapnum_store,
    rx_mapnum,
    aa_mapnum,
):
    """Perform cyclization reaction on peptide with given input and filter for the right product.

    Parameters
    ----------

    df_reactions: pandas.core.frame.DataFrame
        DataFrame storing all available cyclization categories.
    cycl_choice_var: str
        variable notation of cyclization finalchoice
    rx_finalchoice: str
        Contain in case of calc_type == fast the final reaction the user intends to perform. In case of calc_type == slow not defined yet.
    chain_reactant: rdkit.Chem.rdchem.Mol
        Copy of RDkit mol object of the current peptide to perform cyclization reaction.
    additional_reactant: int
        == 0, predefined placeholder for in case an additional reactant will be needed.
    highest_mapnum: int
        In current peptide highest given mapnum.
    mapnum_store: list
        List of mapnums acquired during peptide generation so far.
    rx_mapnum: list
        List of mapnums of the involved amino acids.
    aa_mapnum: list
        List of list with mapnums per amino acid involved in intended cyclization reaction.


    Returns
    ------
    chain_cycl: tuple
        Tuple containing generated reaction products from cyclization reaction filtered for not intended products.
    mapnum_store: list
        Updated list of mapnums acquired during peptide generation.
    additional_mapnum: list
        List with new mapnums for additional component.

    """

    rx_toperform = rdChemReactions.ReactionFromSmarts(
        df_reactions[
            (df_reactions["Type"] == cycl_choice_var)
            & (df_reactions["Variable"] == rx_finalchoice)
        ]
        .reset_index(drop=True)
        .loc[0, "Reaction_Smarts"]
    )
    additional_mapnum = []
    if rx_toperform.GetNumReactantTemplates() == 1:
        chain_cycl = rx_toperform.RunReactant(chain_reactant, additional_reactant)
    elif rx_toperform.GetNumReactantTemplates() == 2:
        additional_reactant = Chem.MolFromSmiles(
            df_reactions[df_reactions["Variable"] == rx_finalchoice]
            .reset_index(drop=True)
            .loc[0, "Second_Reactant"]
        )
        [
            a.SetIntProp("mapnum", highest_mapnum + j + 1)
            for j, a in enumerate(additional_reactant.GetAtoms())
        ]
        highest_mapnum += additional_reactant.GetNumAtoms()
        mapnum_store.append(highest_mapnum)
        rx_mapnum += list(range(mapnum_store[-2] + 1, mapnum_store[-1] + 1))
        additional_mapnum = list(range(mapnum_store[-2] + 1, mapnum_store[-1] + 1))
        chain_cycl = rx_toperform.RunReactants((chain_reactant, additional_reactant))

    chain_cycl = tuple(
        (value,)
        for value in [
            keepmapnum(chain_reactant, additional_reactant, a[0], rx_toperform)
            for a in chain_cycl
        ]
    )
    chain_cycl = find_mapnum_product(chain_cycl, rx_mapnum)
    chain_cycl = remove_intraaacycl_product(chain_cycl, aa_mapnum)

    return chain_cycl, mapnum_store, additional_mapnum


def perform_mod_rx(
    rx_toperform,
    chain_reactant,
    additional_reactant,
    highest_mapnum_mod,
    mapnum_store,
    rx_mapnum,
    aamod_finalchoice,
    df_sequence_overview,
    aa_choice_num,
    famod_mapnum,
):
    """Perform modification reaction on peptide with given input and filter for the right product.

    Parameters
    ----------

    rx_toperform: rdkit.Chem.rdChemReactions.ChemicalReaction
        Modification reaction to perform given by userinput.
    chain_reactant: rdkit.Chem.rdchem.Mol
        Copy of RDkit mol object of the current peptide to perform cyclization reaction.
    additional_reactant: rdkit.Chem.rdchem.Mol
        RDkit mol of additional reactant if needed.
    highest_mapnum_mod: int
        In current peptide highest given mapnum.
    mapnum_store: list
        List of mapnums acquired during peptide generation so far.
    rx_mapnum: list
        List of mapnums of the involved amino acids.
    aamod_finalchoice: str
        Chosen cyclization overcategory (Fatty acid derivatization, Glycosylation, Phosphorylation, Betaalanine, Aspartimide/Isoaspartate, Sulfur Oxidation)
    df_sequence_overview: pd.DataFrame
        DataFrame storing linear sequence information including mapnum.
    aa_choice_num: int
        Number of amino acid to modify.
    famod_mapnum: list
        List containing in case of fatty acid derivatization the mapnum of the fatty acid to attach.

    Returns
    ------
    chain_mod: tuple
        Tuple containing generated reaction products from modification reaction filtered for not intended products.
    mapnum_store: list
        Updated list of mapnums acquired during peptide generation.
    rx_mapnum: list
        Updated list of mapnums of the involved amino acids.

    """

    if rx_toperform.GetNumReactantTemplates() == 1:
        chain_mod = rx_toperform.RunReactant(chain_reactant, additional_reactant)
    elif rx_toperform.GetNumReactantTemplates() == 2:
        chain_mod = rx_toperform.RunReactants((chain_reactant, additional_reactant))
        highest_mapnum = highest_mapnum_mod
        mapnum_store.append(highest_mapnum)
        rx_mapnum = rx_mapnum + list(range(mapnum_store[-2] + 1, mapnum_store[-1] + 1))
    chain_mod = tuple(
        (value,)
        for value in [
            keepmapnum(chain_reactant, additional_reactant, a[0], rx_toperform)
            for a in chain_mod
        ]
    )
    if aamod_finalchoice == "Fatty Acid Derivatization":
        famod_rx_mapnum = df_sequence_overview.loc[aa_choice_num, "Mapnum"] + list(
            range(famod_mapnum[-2] + 1, famod_mapnum[-1] + 1)
        )
        chain_mod = find_mapnum_product(chain_mod, famod_rx_mapnum)
    else:
        chain_mod = find_mapnum_product(chain_mod, rx_mapnum)
    return chain_mod, mapnum_store, rx_mapnum


def assign_rx_aamod_ht_input(aamod_rx_info, df_modifications, errormessage):
    """Assigns for amino acid modification needed variables from the ht input for amino acid modification to perform the intended modification reaction.

    Parameters
    ----------
    aamod_rx_info: list
        List containing number of amino acid to modify and information for the modification reaction.
    df_modifications: pandas.core.frame.Dataframe
        Contains all possible modification reactions.
    errormessage: str
        Current errormessage (none) which can be updated in case an error in reaction assignment for amino acid modification.

    Returns
    ------
    errormessage: str
        Updated errormessage.
    option_choice: str
        String of the reaction option chosen.
    aamod_finalchoice: str
        String of the reaction option chosen. In case of fast calc method this variable corresponds to the option_choice.
    fa_choice: str
        String of fatty acid chosen for modifiation in case of fatty acid modification.
    lnk_choice: str
        String of linker chosen for modifiation in case of fatty acid modification.
    sp_choice: str
        String of spacer chosen for modifiation in case of fatty acid modification.
    sp_num: int
        Number for the times the spacer should be included in case of fatty acid modification. Zero per default.
    lnk_num: int
        Number for the times the linker should be included in case of fatty acid modification. Zero per default.
    """

    option_choice = ""
    aamod_finalchoice = ""
    fa_choice = lnk_choice = sp_choice = ""
    sp_num = lnk_num = 0

    glyc_check = "n_" + aamod_rx_info[1]
    if (
        aamod_rx_info[1] in df_modifications["Variable"].values
        and df_modifications[df_modifications["Variable"] == aamod_rx_info[1]][
            "Type"
        ].iloc[0]
        == "fattyac"
    ):
        if len(aamod_rx_info) == 2:
            aamod_rx_info += [[[1, aamod_rx_info[1]]]]
            aamod_rx_info[1] = "fa"
        elif len(aamod_rx_info) == 3:
            aamod_rx_info[2] += [[1, aamod_rx_info[1]]]
            aamod_rx_info[1] = "fa"
    if (
        aamod_rx_info[1] in df_modifications["Variable"].values
        and df_modifications[df_modifications["Variable"] == aamod_rx_info[1]][
            "ModType"
        ].iloc[0]
        != "fa"
    ):
        aamod_finalchoice = df_modifications[
            df_modifications["Variable"] == aamod_rx_info[1]
        ]["ModType"].iloc[0]
        option_choice = aamod_rx_info[1]
    elif glyc_check in df_modifications["Variable"].values:
        aamod_finalchoice = df_modifications[
            df_modifications["Variable"] == glyc_check
        ]["ModType"].iloc[0]
        aamod_rx_info[1] = glyc_check
        option_choice = aamod_rx_info[1]
    elif aamod_rx_info[1] == "fa":
        aamod_finalchoice = "fa"

        if len(aamod_rx_info[2]) > 3 or len(aamod_rx_info[2]) == 0:
            errormessage = "Fatty acid modification reaction unclear"
        for fa_elem in aamod_rx_info[2]:
            if (
                fa_elem[1]
                in df_modifications[df_modifications["Type"] == "fattyac"][
                    "Variable"
                ].values
                and fa_choice == ""
                and fa_elem[0] == 1
            ):
                fa_choice = df_modifications[
                    df_modifications["Variable"] == fa_elem[1]
                ]["Name"].iloc[0]
            elif (
                fa_elem[1]
                in df_modifications[df_modifications["Type"] == "linker"][
                    "Variable"
                ].values
                and lnk_choice == ""
            ):
                lnk_choice = df_modifications[
                    df_modifications["Variable"] == fa_elem[1]
                ]["Name"].iloc[0]
                lnk_num = fa_elem[0]
            elif (
                fa_elem[1]
                in df_modifications[df_modifications["Type"] == "spacer"][
                    "Variable"
                ].values
                and sp_choice == ""
            ):
                sp_choice = df_modifications[
                    df_modifications["Variable"] == fa_elem[1]
                ]["Name"].iloc[0]
                sp_num = fa_elem[0]
            else:
                errormessage = "Fatty acid information " + fa_elem[1] + " unclear"

        if fa_choice == "":
            errormessage = "Fatty acid information missing"
        if lnk_choice == "":
            lnk_choice = "y-Glu"
        if sp_choice == "":
            sp_choice = "Ado"
    else:
        errormessage = "Modification reaction " + aamod_rx_info[1] + " unclear"

    return (
        errormessage,
        option_choice,
        aamod_finalchoice,
        fa_choice,
        lnk_choice,
        sp_choice,
        sp_num,
        lnk_num,
    )


def fa_assembly(
    df_modifications, fa_choice, lnk_num, lnk_choice, sp_num, sp_choice, highest_mapnum
):
    """Processes user input for fatty acid to attach and assembles the parts of the fatty acid to attach.

    Parameters
    ----------

    df_modifications: pandas.core.frame.Dataframe
        Contains all possible modification reactions.
    fa_choice: str
        String of fatty acid chosen for modifiation in case of fatty acid modification.
    lnk_num: int
        Number for the times the linker should be included in case of fatty acid modification.
    lnk_choice: str
        String of linker chosen for modifiation in case of fatty acid modification.
    sp_num: int
        Number for the times the spacer should be included in case of fatty acid modification.
    sp_choice: str
        String of spacer chosen for modifiation in case of fatty acid modification.
    highest_mapnum: int
        Currently highest mapnum in the peptide.


    Returns
    ------
    famod: rdkit.Chem.rdchem.Mol
        RDkit Mol of the fatty acid to attach to the peptide.
    fa_info: str
        String containing the information of the fatty acid to attach.
    famod_highest_mapnum: int
        Highest mapnum in the fatty acid.
    famod_mapnum: list
        List of mapnums added to the fatty acid to attach continuing from the highest mapnum in the peptide.

    """

    if int(lnk_num) == 0:
        lnk_info = ""
    elif int(lnk_num) == 1:
        lnk_info = lnk_choice + "-"
    else:
        lnk_info = str(lnk_num) + " " + lnk_choice + "-"

    if int(sp_num) == 0:
        sp_info = ""
    elif int(sp_num) == 1:
        sp_info = sp_choice + "-"
    else:
        sp_info = str(sp_num) + " " + sp_choice + "-"

    fa_info = "*" + sp_info + lnk_info + fa_choice

    fa_smiles = [
        df_modifications[df_modifications["Name"] == fa_choice].iloc[0]["Smiles"]
    ]
    lnk_smiles = [
        df_modifications[df_modifications["Name"] == lnk_choice].iloc[0]["Smiles"]
    ]
    lnk_rxn = [
        df_modifications[df_modifications["Name"] == lnk_choice].iloc[0][
            "SmartsReaction"
        ]
    ]
    sp_smiles = [
        df_modifications[df_modifications["Name"] == sp_choice].iloc[0]["Smiles"]
    ]
    sp_rxn = [
        df_modifications[df_modifications["Name"] == sp_choice].iloc[0][
            "SmartsReaction"
        ]
    ]

    lst_famod_smiles = fa_smiles + lnk_smiles * int(lnk_num) + sp_smiles * int(sp_num)
    lst_famod_rx = lnk_rxn * int(lnk_num) + sp_rxn * int(sp_num)
    famod_mols = [Chem.MolFromSmiles(i) for i in lst_famod_smiles]

    famod = famod_mols[0]
    famod_highest_mapnum = highest_mapnum
    famod_mapnum = [famod_highest_mapnum]
    [
        a.SetIntProp("mapnum", famod_highest_mapnum + j + 1)
        for j, a in enumerate(famod.GetAtoms())
    ]
    famod_highest_mapnum += famod.GetNumAtoms()
    famod_mapnum.append(famod_highest_mapnum)

    for i in range(1, len(famod_mols)):
        atc = famod_mols[i]
        [a.ClearProp("old_mapno") for a in famod.GetAtoms()]
        [
            a.SetIntProp("mapnum", famod_highest_mapnum + j + 1)
            for j, a in enumerate(atc.GetAtoms())
        ]
        famod_highest_mapnum += atc.GetNumAtoms()
        famod_mapnum.append(famod_highest_mapnum)
        famod_current_mapnum = list(
            range(famod_mapnum[i - 1] + 1, famod_mapnum[i + 1] + 1)
        )
        faassembly_rx_smarts = lst_famod_rx[i - 1]
        faassembly_rx = rdChemReactions.ReactionFromSmarts(faassembly_rx_smarts)
        fa_chain_mod = faassembly_rx.RunReactants((famod, atc))

        fa_chain_mod = tuple(
            (value,)
            for value in [
                keepmapnum(famod, atc, a[0], faassembly_rx) for a in fa_chain_mod
            ]
        )
        fa_chain_mod = find_mapnum_product(fa_chain_mod, famod_current_mapnum)

        fa_chain_ext = fa_chain_mod[0][0]
        Chem.SanitizeMol(fa_chain_ext)
        famod = fa_chain_ext

    return famod, fa_info, famod_highest_mapnum, famod_mapnum


# %%
## IMPORT AMINOACIDS/N-&C-TERMINUS/CYCLIZATION REACTIONS/MODIFICATION REACTIONS ##

if __name__ == "__main__":

    df_aminoacids = (pd.read_csv(path_to_AminoAcids)).fillna(0)

    df_termmod = (pd.read_csv(path_to_TerminalMod)).fillna(0)
    dict_nterm = (
        df_termmod[df_termmod["Mod_Site"] == "N"]
        .set_index("Entry_Code")
        .to_dict("index")
    )
    dict_cterm = (
        df_termmod[df_termmod["Mod_Site"] == "C"]
        .set_index("Entry_Code")
        .to_dict("index")
    )

    df_reactions = (pd.read_csv(path_to_CyclReactions)).fillna(0)
    dict_rxtypes = {
        "headtail": "Head to Tail",
        "headsc": "Head to Sidechain",
        "tailsc": "Tail to Sidechain",
        "scsc": "Sidechain to Sidechain",
        "twocomp": "Two Component",
        "bicycle": "Bicycle Formation",
    }
    df_modifications = (pd.read_csv(path_to_ModReactions)).fillna(0)

    with open(path_to_dict_default, "rb") as file:
        dict_default = pickle.load(file)

    try:

        nterm_opt = ["none"] + list(dict_nterm.keys())
        cterm_opt = ["none"] + list(dict_cterm.keys())
        allowed = r"[^A-Za-z0-9]"

        window_status = 1
        back = cycl = aamod = False
        default_sequence = " PICKAPEP "

        while window_status == 1:
            # Linear peptide assembly
            (
                sequence,
                cterm,
                nterm,
                cycl,
                window_status,
                default_pept,
                df_aminoacids,
                calc_type,
                calc_table,
                calc_table_path,
            ) = userinput1(
                nterm_opt,
                cterm_opt,
                dict_default,
                window_status,
                df_aminoacids,
                default_sequence,
                path_to_teriyaki,
            )

            if window_status == 4:
                break

            skip = 0
            fastcalc_pepnum = 0
            errormessage = "none"

            while (fastcalc_pepnum == 0 or fastcalc_pepnum <= len(calc_table)) and (
                window_status != 1 and window_status != 4
            ):

                pepdone = False
                fastcalc_pepnum += 1
                if calc_type == "fast":
                    cycl = 1
                    if fastcalc_pepnum - 1 == len(calc_table):
                        window_status = 3
                        curr_status = calc_table.loc[fastcalc_pepnum - 2, "Status"]
                        if curr_status == "Done" or curr_status == errormessage:
                            if curr_status == "Done":
                                skip += 1
                        pepdone = True

                if default_pept == "" and pepdone == False:

                    errormessage = "none"

                    if calc_type == "fast":
                        sequence = calc_table.at[
                            fastcalc_pepnum - 1, "Input Linear Sequence"
                        ]

                    aaseq, aa_smiles, errormessage = check_sequence_input(
                        sequence, allowed, df_aminoacids
                    )
                    if calc_type == "slow":
                        default_sequence = "".join(aaseq)

                    if aaseq == []:
                        errormessage = ""

                    if errormessage != "none":
                        window_status, calc_table = errorupdate(
                            window_status,
                            calc_type,
                            calc_table,
                            fastcalc_pepnum,
                            errormessage,
                            new_window_status=1,
                        )
                        continue

                    else:
                        df_sequence_overview = pd.DataFrame(
                            columns=["AminoAcid", "Mapnum", "AminoAcidInfo"]
                        )
                        df_cycl_overview = pd.DataFrame(
                            columns=[
                                "AA1",
                                "Num AA1",
                                "AA2",
                                "Num AA2",
                                "AA3",
                                "NumAA3",
                                "CyclMyMap",
                                "CyclType",
                                "CyclName",
                                "Additional",
                            ]
                        )
                        df_modif_overview = pd.DataFrame(
                            columns=[
                                "AAModif",
                                "NumAAModif",
                                "ModifMyMap",
                                "ModifType",
                                "Additional",
                            ]
                        )

                        (
                            chain,
                            df_sequence_overview,
                            mapnum_store,
                            highest_mapnum,
                            errormessage,
                        ) = create_linear_peptide(
                            aaseq, aa_smiles, df_sequence_overview, errormessage
                        )

                        if errormessage != "none":
                            window_status, calc_table = errorupdate(
                                window_status,
                                calc_type,
                                calc_table,
                                fastcalc_pepnum,
                                errormessage,
                                new_window_status=1,
                            )
                            continue
                        else:
                            # get nterm and cterm input from table in highthroughput mode
                            if calc_type == "fast":
                                nterm = calc_table.at[
                                    fastcalc_pepnum - 1, "Input Nterm"
                                ]
                                nterm = re.sub(r"[^A-Za-z]", "", str(nterm))
                                if (
                                    nterm.lower() == "none"
                                    or nterm == " "
                                    or nterm == ""
                                ):
                                    nterm = "none"
                                else:
                                    nterm = nterm[0].upper() + nterm[1:].lower() + "-"
                                if nterm not in nterm_opt:
                                    errormessage = "N-Terminus not available"
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=1,
                                    )
                                    continue
                                cterm = calc_table.at[
                                    fastcalc_pepnum - 1, "Input Cterm"
                                ]
                                cterm = re.sub(r"[^A-Za-z]", "", str(cterm))
                                if (
                                    cterm.lower() == "none"
                                    or cterm == " "
                                    or cterm == ""
                                ):
                                    cterm = "none"
                                else:
                                    cterm = "-" + cterm[0].upper() + cterm[1:].lower()
                                if cterm not in cterm_opt and errormessage == "none":
                                    errormessage = "C-Terminus not available"
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=1,
                                    )
                                    continue

                            ## Modify N-terminus
                            if nterm != "none":
                                (
                                    chain,
                                    df_sequence_overview,
                                    mapnum_store,
                                    highest_mapnum,
                                    errormessage,
                                ) = attach_nterminus(
                                    chain,
                                    df_sequence_overview,
                                    mapnum_store,
                                    highest_mapnum,
                                    dict_nterm,
                                    nterm,
                                    errormessage,
                                )
                                if errormessage != "none":
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=1,
                                    )
                                    continue
                            else:
                                df_sequence_overview.at[0, "Mapnum"] = (
                                    df_sequence_overview.at[1, "Mapnum"]
                                )
                                df_sequence_overview.at[0, "AminoAcid"] = "none"

                            ## Modify C-terminus
                            if cterm != "none":
                                (
                                    chain,
                                    df_sequence_overview,
                                    mapnum_store,
                                    highest_mapnum,
                                    errormessage,
                                ) = attach_cterminus(
                                    chain,
                                    df_sequence_overview,
                                    mapnum_store,
                                    highest_mapnum,
                                    dict_cterm,
                                    cterm,
                                    errormessage,
                                )
                                if errormessage != "none":
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=1,
                                    )
                                    continue
                            else:
                                df_sequence_overview.at[
                                    len(df_sequence_overview) - 1, "Mapnum"
                                ] = df_sequence_overview.at[
                                    len(df_sequence_overview) - 2, "Mapnum"
                                ]
                                df_sequence_overview.at[
                                    len(df_sequence_overview) - 1, "AminoAcid"
                                ] = "none"

                            [a.ClearProp("old_mapno") for a in chain.GetAtoms()]

                            # Extend terminus on one atom to ensure valid reaction template for head/tail involved cyclization reactions
                            nterm_atoms_idx = [
                                atom.GetIdx()
                                for atom in chain.GetAtoms()
                                if atom.HasProp("mapnum")
                                and atom.GetIntProp("mapnum")
                                in df_sequence_overview["Mapnum"].iloc[0]
                            ]
                            cterm_atoms_idx = [
                                atom.GetIdx()
                                for atom in chain.GetAtoms()
                                if atom.HasProp("mapnum")
                                and atom.GetIntProp("mapnum")
                                in df_sequence_overview["Mapnum"].iloc[-1]
                            ]

                            for idx in nterm_atoms_idx:
                                atom_neighbors = (
                                    chain.GetAtomWithIdx(idx)
                                ).GetNeighbors()
                                new_neighbor_idx = [
                                    i.GetIdx()
                                    for i in atom_neighbors
                                    if i.GetIdx() not in nterm_atoms_idx
                                ]
                                df_sequence_overview["Mapnum"].iloc[
                                    0
                                ] = df_sequence_overview["Mapnum"].iloc[0] + [
                                    chain.GetAtomWithIdx(i).GetIntProp("mapnum")
                                    for i in new_neighbor_idx
                                ]

                            for idx in cterm_atoms_idx:
                                atom_neighbors = (
                                    chain.GetAtomWithIdx(idx)
                                ).GetNeighbors()
                                new_neighbor_idx = [
                                    i.GetIdx()
                                    for i in atom_neighbors
                                    if i.GetIdx() not in cterm_atoms_idx
                                ]
                                df_sequence_overview["Mapnum"].iloc[
                                    -1
                                ] = df_sequence_overview["Mapnum"].iloc[-1] + [
                                    chain.GetAtomWithIdx(i).GetIntProp("mapnum")
                                    for i in new_neighbor_idx
                                ]

                            MyLinearPeptide = Peptide(
                                chain,
                                mapnum_store,
                                df_sequence_overview,
                                df_cycl_overview,
                                df_modif_overview,
                                assign_pepinfo(chain, df_sequence_overview),
                            )
                            peptide_memory = [MyLinearPeptide]

                else:
                    if calc_type == "fast":
                        pass
                    else:
                        peptide_memory = dict_default[default_pept].copy()
                        # get parameters from default peptide (class) to allow futher processing
                        current_peptide = peptide_memory[-1]
                        chain = current_peptide.peptide
                        highest_mapnum = current_peptide.highestmapnum
                        mapnum_store = current_peptide.mapnumstore.copy()
                        df_sequence_overview = current_peptide.df_sequence.copy()
                        df_cycl_overview = current_peptide.df_cycl.copy()
                        df_modif_overview = current_peptide.df_aamod.copy()
                        cycl_count = current_peptide.cycl_num
                        aamod_count = current_peptide.aamod_num
                        df_pepinfo = current_peptide.pepinfo.copy()
                        aa_count = len(df_sequence_overview) - 2

                cyclization_choices = list(dict_rxtypes.values())
                modification_choices = [
                    "Fatty Acid Derivatization",
                    "Glycosylation",
                    "Phosphorylation",
                    "Betaalanine Formation",
                    "Aspartimide / Isoaspartate",
                    "Sulfur Oxidation",
                ]  # change as dict calc above
                cycl_count = 0
                aamod_count = 0
                cycl_ht_todo = 0
                aamod_ht_todo = 0

                while window_status == 2 or window_status == 3:
                    if window_status == 2:

                        if cycl or aamod or back:

                            current_peptide = peptide_memory[-1]
                            chain = current_peptide.peptide
                            highest_mapnum = current_peptide.highestmapnum
                            mapnum_store = current_peptide.mapnumstore.copy()
                            df_sequence_overview = current_peptide.df_sequence.copy()
                            df_cycl_overview = current_peptide.df_cycl.copy()
                            df_modif_overview = current_peptide.df_aamod.copy()
                            cycl_count = current_peptide.cycl_num
                            aamod_count = current_peptide.aamod_num
                            aa_count = (
                                len(df_sequence_overview) - 2
                            )  # correct for N-&C-terminus in table

                            # Sort cycl/aamod input either in high throughput mode or get userinput for cyclization/aamod
                            if calc_type == "fast":
                                cycl_htinput = calc_table.at[
                                    fastcalc_pepnum - 1, "Input Cyclizations"
                                ]
                                cycl_htinput_sorted = sort_cycl_htinput(
                                    cycl_htinput, df_sequence_overview
                                )
                                cycl_ht_todo = len(cycl_htinput_sorted)
                                aamod_htinput = calc_table.at[
                                    fastcalc_pepnum - 1, "Input Modifications"
                                ]
                                aamod_htinput_sorted = sort_aamod_htinput(
                                    aamod_htinput,
                                    df_sequence_overview,
                                    df_modifications,
                                )
                                aamod_ht_todo = len(aamod_htinput_sorted)
                                if cycl_ht_todo == 0 and aamod_ht_todo == 0:
                                    cycl = False
                                    aamod = False
                                else:
                                    if cycl_ht_todo > 0:
                                        cycl = True
                                    if aamod_ht_todo > 0:
                                        aamod = True

                            else:
                                errormessage = "none"
                                (
                                    cycl,
                                    cyclization_finalchoice,
                                    aamod,
                                    aamod_finalchoice,
                                    peptide_memory,
                                    back,
                                    window_status,
                                ) = userinput2(
                                    cyclization_choices,
                                    modification_choices,
                                    df_sequence_overview["AminoAcidInfo"].tolist(),
                                    peptide_memory,
                                    window_status,
                                    path_to_teriyaki,
                                    path_to_teriyaki_easteregg,
                                )
                                if back != True:
                                    cycl_ht_todo = cycl_count + 1
                                    aamod_ht_todo = aamod_count + 1

                            while cycl_count < cycl_ht_todo and cycl:

                                current_peptide = peptide_memory[-1]
                                chain = current_peptide.peptide
                                highest_mapnum = current_peptide.highestmapnum
                                mapnum_store = current_peptide.mapnumstore.copy()
                                df_sequence_overview = (
                                    current_peptide.df_sequence.copy()
                                )
                                df_cycl_overview = current_peptide.df_cycl.copy()
                                df_modif_overview = current_peptide.df_aamod.copy()
                                cycl_count = current_peptide.cycl_num
                                aamod_count = current_peptide.aamod_num
                                aa_count = (
                                    len(df_sequence_overview) - 2
                                )  # correct for N-&C-terminus in table

                                rx_possible = False
                                [a.ClearProp("old_mapno") for a in chain.GetAtoms()]

                                if calc_type == "fast":

                                    additional_needed = False
                                    additional_output = ""
                                    cycl_rx_info = cycl_htinput_sorted[cycl_count]

                                    if len(cycl_rx_info[0]) == 2:
                                        if any(
                                            num == 0 for num in cycl_rx_info[0]
                                        ) and any(
                                            num == len(df_sequence_overview) - 1
                                            for num in cycl_rx_info[0]
                                        ):
                                            cycl_choice_var = "headtail"
                                            cycl_rx_info[1] = "amide"
                                        elif any(num == 0 for num in cycl_rx_info[0]):
                                            cycl_choice_var = "tailsc"
                                            cycl_rx_info[1] = "amide"
                                        elif any(
                                            num == len(df_sequence_overview) - 1
                                            for num in cycl_rx_info[0]
                                        ):
                                            cycl_choice_var = "headsc"
                                            cycl_rx_info[1] = "amide"
                                        elif df_reactions.index[
                                            (
                                                df_reactions["Variable"].str.lower()
                                                == cycl_rx_info[1].lower()
                                            )
                                            & (df_reactions["Type"] == "twocomp")
                                        ].tolist():
                                            cycl_choice_var = "twocomp"
                                        else:
                                            cycl_choice_var = "scsc"
                                    elif len(cycl_rx_info[0]) == 3:
                                        cycl_choice_var = "bicycle"
                                    else:
                                        errormessage = (
                                            "Reaction "
                                            + str(cycl_count + 1)
                                            + " unclear aminoacid input"
                                        )
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        break
                                    matching_index = df_reactions.index[
                                        (
                                            df_reactions["Variable"].str.lower()
                                            == cycl_rx_info[1].lower()
                                        )
                                        & (df_reactions["Type"] == cycl_choice_var)
                                    ].tolist()

                                    if len(matching_index) == 0:
                                        errormessage = (
                                            "Reaction "
                                            + cycl_rx_info[1]
                                            + " not found!"
                                        )
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        break
                                    else:
                                        cyclization_finalchoice = dict_rxtypes[
                                            cycl_choice_var
                                        ]
                                        aa_choice_num_output = cycl_rx_info[0]
                                        aa_choice_output = []
                                        aa_options = df_sequence_overview.loc[
                                            :, "AminoAcidInfo"
                                        ].tolist()
                                        for i in aa_choice_num_output:
                                            aa_info = aa_options[i]
                                            aa_choice_output.append(aa_info)
                                        rx_finalchoice = df_reactions.loc[
                                            matching_index[0], "Variable"
                                        ]
                                else:
                                    (
                                        aa_choice_output,
                                        aa_choice_num_output,
                                        additional_needed,
                                        additional_info,
                                        additional_output,
                                        skip_cycl,
                                    ) = sort_cycl_input(
                                        df_sequence_overview,
                                        aa_count,
                                        cyclization_finalchoice,
                                        df_reactions,
                                        userinput2_cycl,
                                    )
                                    cycl_choice_var = next(
                                        (
                                            reaction
                                            for reaction, value in dict_rxtypes.items()
                                            if value == cyclization_finalchoice
                                        ),
                                        None,
                                    )
                                    rx_finalchoice = ""
                                    if skip_cycl:
                                        back = True
                                        cycl_ht_todo -= 1
                                        break

                                chain_reactant = chain
                                additional_reactant = 0

                                possible_rx, rx_mapnum, aa_mapnum, additional_output = (
                                    get_possible_cycl_rx(
                                        aa_choice_num_output,
                                        df_sequence_overview,
                                        df_reactions,
                                        chain,
                                        cyclization_finalchoice,
                                        rx_finalchoice,
                                        cycl_choice_var,
                                        additional_output,
                                        calc_type,
                                    )
                                )

                                if calc_type == "fast":
                                    if rx_finalchoice in possible_rx:
                                        rx_possible = True
                                    else:
                                        errormessage = (
                                            "Reaction "
                                            + str(cycl_count + 1)
                                            + " not possible"
                                        )
                                else:
                                    if len(possible_rx) > 1:
                                        possible_rx_name = [
                                            df_reactions.loc[
                                                df_reactions["Variable"] == rx_var,
                                                "Name",
                                            ].iloc[0]
                                            for rx_var in possible_rx
                                        ]
                                        rx_finalchoice_name = input_singlechoice(
                                            possible_rx_name
                                        )
                                        if rx_finalchoice_name == "":
                                            back = True
                                            cycl_ht_todo -= 1
                                            break
                                        else:
                                            rx_finalchoice = df_reactions.loc[
                                                df_reactions["Name"]
                                                == rx_finalchoice_name,
                                                "Variable",
                                            ].iloc[0]
                                            rx_possible = True
                                    elif len(possible_rx) == 1:
                                        rx_finalchoice = possible_rx[0]
                                        rx_possible = True
                                    else:
                                        errormessage = "Reaction not possible"

                                if rx_possible:

                                    chain_cycl, mapnum_store, additional_mapnum = (
                                        perform_cycl_rx(
                                            df_reactions,
                                            cycl_choice_var,
                                            rx_finalchoice,
                                            chain_reactant,
                                            additional_reactant,
                                            highest_mapnum,
                                            mapnum_store,
                                            rx_mapnum,
                                            aa_mapnum,
                                        )
                                    )

                                    if num_equiv_solutions(chain_cycl) == 1:
                                        n = 0
                                    else:
                                        if calc_type == "fast":
                                            errormessage = (
                                                "Reaction "
                                                + str(cycl_count + 1)
                                                + " results in "
                                                + str(num_equiv_solutions(chain_cycl))
                                                + " products"
                                            )
                                        else:
                                            solution_ignore = tk.messagebox.askyesno(
                                                "Question",
                                                "The reaction would result in "
                                                + str(num_equiv_solutions(chain_cycl))
                                                + " products. Do you want to proceed with one of the products?",
                                            )
                                            if solution_ignore:
                                                n = random.randrange(
                                                    0,
                                                    num_equiv_solutions(chain_cycl),
                                                    1,
                                                )
                                            else:
                                                back = True
                                                cycl_ht_todo -= 1
                                                break

                                    if errormessage == "none":
                                        chain_ext = chain_cycl[n][0]
                                        Chem.SanitizeMol(chain_ext)
                                        chain = chain_ext
                                        rx_name = df_reactions.loc[
                                            df_reactions["Variable"] == rx_finalchoice,
                                            "Name",
                                        ].iloc[0]
                                        df_cycl_overview.loc[
                                            cycl_count,
                                            [
                                                "AA1",
                                                "Num AA1",
                                                "AA2",
                                                "Num AA2",
                                                "CyclMyMap",
                                                "CyclType",
                                                "CyclName",
                                                "AA3",
                                                "NumAA3",
                                                "Additional",
                                            ],
                                        ] = [
                                            aa_choice_output[0],
                                            aa_choice_num_output[0],
                                            aa_choice_output[1],
                                            aa_choice_num_output[1],
                                            rx_mapnum,
                                            cyclization_finalchoice,
                                            rx_name,
                                            0,
                                            0,
                                            0,
                                        ]
                                        if additional_needed:
                                            df_cycl_overview.loc[
                                                cycl_count, ["Additional"]
                                            ] = (additional_info + additional_output)
                                        if (
                                            cyclization_finalchoice
                                            == "Bicycle Formation"
                                        ):
                                            df_cycl_overview.loc[
                                                cycl_count, ["AA3", "NumAA3"]
                                            ] = [
                                                aa_choice_output[2],
                                                aa_choice_num_output[2],
                                            ]
                                        for i in aa_choice_num_output:
                                            df_sequence_overview.loc[
                                                i, "AminoAcidInfo"
                                            ] += (
                                                " * cycl"
                                                + "("
                                                + str(cycl_count + 1)
                                                + ")"
                                            )
                                            df_sequence_overview.at[i, "Mapnum"] = list(
                                                set(
                                                    df_sequence_overview.at[i, "Mapnum"]
                                                    + additional_mapnum
                                                )
                                            )
                                        MyModPeptide = Peptide(
                                            chain,
                                            mapnum_store,
                                            df_sequence_overview,
                                            df_cycl_overview,
                                            df_modif_overview,
                                            assign_pepinfo(chain, df_sequence_overview),
                                        )
                                        peptide_memory.append(MyModPeptide)
                                        cycl_count += 1
                                    else:
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        back = True
                                        cycl = False
                                        if fastcalc_pepnum == len(calc_table):
                                            window_status = 3
                                        break
                                else:
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=2,
                                    )
                                    cycl = False
                                    back = True
                                    if fastcalc_pepnum == len(calc_table):
                                        window_status = 3
                                    break

                            while aamod_count < aamod_ht_todo and aamod:

                                current_peptide = peptide_memory[-1]
                                chain = current_peptide.peptide
                                highest_mapnum = current_peptide.highestmapnum
                                mapnum_store = current_peptide.mapnumstore.copy()
                                df_sequence_overview = (
                                    current_peptide.df_sequence.copy()
                                )
                                df_cycl_overview = current_peptide.df_cycl.copy()
                                df_modif_overview = current_peptide.df_aamod.copy()
                                cycl_count = current_peptide.cycl_num
                                aamod_count = current_peptide.aamod_num
                                aa_count = (
                                    len(df_sequence_overview) - 2
                                )  # correct for N-&C-terminus in table

                                rx_possible = False
                                [a.ClearProp("old_mapno") for a in chain.GetAtoms()]

                                highest_mapnum_mod = famod_mapnum = 0

                                aa_options = df_sequence_overview.loc[
                                    1:aa_count, "AminoAcidInfo"
                                ].tolist()

                                if calc_type == "fast":
                                    aamod_rx_info = aamod_htinput_sorted[aamod_count]
                                    if aamod_rx_info[0] == []:
                                        errormessage = (
                                            "Aminoacid position for modification "
                                            + str(aamod_count + 1)
                                            + " unclear"
                                        )
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        break
                                    elif (
                                        aamod_rx_info[0][0] >= 1
                                        and aamod_rx_info[0][0] <= aa_count
                                    ):
                                        aa_choice_num = aamod_rx_info[0][0]
                                        aa_choice = aa_options[aa_choice_num - 1]
                                    else:
                                        errormessage = (
                                            "Aminoacid position for modification "
                                            + str(aamod_count + 1)
                                            + " unclear"
                                        )
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        break

                                    (
                                        errormessage,
                                        option_choice,
                                        aamod_finalchoice,
                                        fa_choice,
                                        lnk_choice,
                                        sp_choice,
                                        sp_num,
                                        lnk_num,
                                    ) = assign_rx_aamod_ht_input(
                                        aamod_rx_info, df_modifications, errormessage
                                    )

                                    option_choice_var = ""

                                    if errormessage != "none":
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        back = True
                                        aamod = False
                                        if fastcalc_pepnum == len(calc_table):
                                            window_status = 3
                                        break

                                    if option_choice:
                                        option_choice_var = option_choice
                                        option_choice = (
                                            df_modifications[
                                                df_modifications["Variable"]
                                                == option_choice_var
                                            ]
                                            .reset_index(drop=True)
                                            .loc[0, "Name"]
                                        )

                                chain_reactant = chain
                                additional_reactant = 0

                                # Twocomponent reactions
                                if (
                                    aamod_finalchoice == "Fatty Acid Derivatization"
                                    or aamod_finalchoice == "fa"
                                    or aamod_finalchoice == "Glycosylation"
                                    or aamod_finalchoice == "glyc"
                                ):
                                    if (
                                        aamod_finalchoice == "Fatty Acid Derivatization"
                                        or aamod_finalchoice == "fa"
                                    ):
                                        aamod_finalchoice = "Fatty Acid Derivatization"
                                        fa_options = df_modifications[
                                            df_modifications["Type"] == "fattyac"
                                        ]["Name"].tolist()
                                        lnk_options = df_modifications[
                                            df_modifications["Type"] == "linker"
                                        ]["Name"].tolist()
                                        sp_options = df_modifications[
                                            df_modifications["Type"] == "spacer"
                                        ]["Name"].tolist()

                                        if calc_type != "fast":
                                            (
                                                aa_choice,
                                                fa_choice,
                                                lnk_num,
                                                lnk_choice,
                                                sp_num,
                                                sp_choice,
                                            ) = userinput2_famod(
                                                aa_options,
                                                fa_options,
                                                lnk_options,
                                                sp_options,
                                                aamod_finalchoice,
                                                path_to_teriyaki,
                                            )
                                            if aa_choice == "":
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                            aa_choice_num = int(
                                                [
                                                    part
                                                    for part in aa_choice.split(" ")
                                                    if part
                                                ][0]
                                            )

                                        rx_mapnum = df_sequence_overview.loc[
                                            aa_choice_num, "Mapnum"
                                        ]

                                        (
                                            famod,
                                            fa_info,
                                            famod_highest_mapnum,
                                            famod_mapnum,
                                        ) = fa_assembly(
                                            df_modifications,
                                            fa_choice,
                                            lnk_num,
                                            lnk_choice,
                                            sp_num,
                                            sp_choice,
                                            highest_mapnum,
                                        )

                                        df_fa_mod_smarts = df_modifications[
                                            df_modifications["Type"] == "fattyac_atc"
                                        ]
                                        possible_rx = []
                                        rx_mapnum_atomidx = [
                                            atom.GetIdx()
                                            for atom in chain.GetAtoms()
                                            if atom.GetIntProp("mapnum") in rx_mapnum
                                        ]
                                        for rx in df_fa_mod_smarts["SmartsReaction"]:
                                            x = rdChemReactions.ReactionFromSmarts(rx)
                                            for match in chain.GetSubstructMatches(
                                                x.GetReactantTemplate(0)
                                            ):
                                                if all(
                                                    idx in rx_mapnum_atomidx
                                                    for idx in match
                                                ):
                                                    possible_rx = possible_rx + [x]
                                        if len(possible_rx) >= 1:
                                            rx_toperform = possible_rx[0]
                                        else:
                                            errormessage = (
                                                "Fatty acid attachement not possible"
                                            )
                                            if calc_type == "fast":
                                                window_status, calc_table = errorupdate(
                                                    window_status,
                                                    calc_type,
                                                    calc_table,
                                                    fastcalc_pepnum,
                                                    errormessage,
                                                    new_window_status=2,
                                                )
                                                aamod = False
                                                if fastcalc_pepnum == len(calc_table):
                                                    window_status = 3
                                            else:
                                                messagebox.showwarning(
                                                    "Info", errormessage
                                                )
                                                aamod_ht_todo -= 1
                                            back = True
                                            break
                                        additional_reactant = famod
                                        highest_mapnum_mod = famod_highest_mapnum
                                        [
                                            a.ClearProp("old_mapno")
                                            for a in famod.GetAtoms()
                                        ]

                                    elif (
                                        aamod_finalchoice == "Glycosylation"
                                        or aamod_finalchoice == "glyc"
                                    ):
                                        aamod_finalchoice = "Glycosylation"
                                        if calc_type != "fast":
                                            option_choices_all = df_modifications[
                                                (df_modifications["ModType"] == "glyc")
                                                & (df_modifications["Type"] == "n")
                                            ]["Name"].tolist()
                                            aa_choice, option_choice = (
                                                userinput2_aamod_option(
                                                    aa_options,
                                                    option_choices_all,
                                                    aamod_finalchoice,
                                                    path_to_teriyaki,
                                                )
                                            )
                                            if aa_choice == "":
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                            aa_choice_num = int(
                                                [
                                                    part
                                                    for part in aa_choice.split(" ")
                                                    if part
                                                ][0]
                                            )
                                        rx_mapnum = df_sequence_overview.loc[
                                            aa_choice_num, "Mapnum"
                                        ]
                                        additional_reactant = Chem.MolFromSmiles(
                                            str(
                                                df_modifications[
                                                    (
                                                        df_modifications["Name"]
                                                        == option_choice
                                                    )
                                                    & (df_modifications["Type"] == "n")
                                                ].iloc[0]["Smiles"]
                                            )
                                        )
                                        [
                                            a.SetIntProp(
                                                "mapnum", highest_mapnum + j + 1
                                            )
                                            for j, a in enumerate(
                                                additional_reactant.GetAtoms()
                                            )
                                        ]
                                        highest_mapnum_mod = (
                                            highest_mapnum
                                            + additional_reactant.GetNumAtoms()
                                        )
                                        rx_toperform = (
                                            rdChemReactions.ReactionFromSmarts(
                                                str(
                                                    df_modifications[
                                                        df_modifications["Name"]
                                                        == option_choice
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                        )
                                        rx_mapnum_atomidx = [
                                            atom.GetIdx()
                                            for atom in chain.GetAtoms()
                                            if atom.GetIntProp("mapnum") in rx_mapnum
                                        ]
                                        rx_o_glyc = rdChemReactions.ReactionFromSmarts(
                                            str(
                                                df_modifications[
                                                    (
                                                        df_modifications["Name"]
                                                        == option_choice
                                                    )
                                                    & (
                                                        df_modifications["ModType"]
                                                        == "glyc"
                                                    )
                                                    & (df_modifications["Type"] == "o")
                                                ].iloc[0]["SmartsReaction"]
                                            )
                                        )
                                        for match in chain.GetSubstructMatches(
                                            rx_o_glyc.GetReactantTemplate(0)
                                        ):
                                            if all(
                                                idx in rx_mapnum_atomidx
                                                for idx in match
                                            ):
                                                rx_toperform = rx_o_glyc

                                else:  # one component reactions
                                    if (
                                        aamod_finalchoice == "Sulfur Oxidation"
                                        or aamod_finalchoice == "so"
                                    ):
                                        aamod_finalchoice = "Sulfur Oxidation"
                                        if calc_type != "fast":
                                            option_choices = df_modifications[
                                                df_modifications["ModType"] == "so"
                                            ]["Name"].tolist()
                                            aa_choice, option_choice = (
                                                userinput2_aamod_option(
                                                    aa_options,
                                                    option_choices,
                                                    aamod_finalchoice,
                                                    path_to_teriyaki,
                                                )
                                            )
                                            if aa_choice == "":
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                            aa_choice_num = int(
                                                [
                                                    part
                                                    for part in aa_choice.split(" ")
                                                    if part
                                                ][0]
                                            )
                                            option_choice_var = (
                                                df_modifications[
                                                    df_modifications["Name"]
                                                    == option_choice
                                                ]
                                                .reset_index(drop=True)
                                                .loc[0, "Variable"]
                                            )
                                        rx_toperform = (
                                            rdChemReactions.ReactionFromSmarts(
                                                str(
                                                    df_modifications[
                                                        df_modifications["Variable"]
                                                        == option_choice_var
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                        )
                                        rx_mapnum = df_sequence_overview.loc[
                                            aa_choice_num, "Mapnum"
                                        ]
                                    elif (
                                        aamod_finalchoice
                                        == "Aspartimide / Isoaspartate"
                                        or aamod_finalchoice == "asp"
                                        or aamod_finalchoice == "isoasp"
                                    ):
                                        if calc_type != "fast":
                                            option_choices = [
                                                "Aspartimide Formation",
                                                "Isoaspartate Formation",
                                            ]
                                            aa_choice, option_choice = (
                                                userinput2_aamod_option(
                                                    aa_options,
                                                    option_choices,
                                                    aamod_finalchoice,
                                                    path_to_teriyaki,
                                                )
                                            )
                                            if aa_choice == "":
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                            aa_choice_num = int(
                                                [
                                                    part
                                                    for part in aa_choice.split(" ")
                                                    if part
                                                ][0]
                                            )
                                        if (
                                            option_choice == "Aspartimide Formation"
                                            or aamod_finalchoice == "asp"
                                        ):
                                            option_choice = "Aspartimide Formation"
                                            rx_toperform = (
                                                rdChemReactions.ReactionFromSmarts(
                                                    df_modifications[
                                                        df_modifications["ModType"]
                                                        == "asp"
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                            rx_mapnum = (
                                                df_sequence_overview.loc[
                                                    aa_choice_num, "Mapnum"
                                                ]
                                                + df_sequence_overview.loc[
                                                    aa_choice_num + 1, "Mapnum"
                                                ]
                                            )
                                        elif (
                                            option_choice == "Isoaspartate Formation"
                                            or aamod_finalchoice == "isoasp"
                                        ):
                                            option_choice = "Isoaspartate Formation"
                                            rx_toperform = (
                                                rdChemReactions.ReactionFromSmarts(
                                                    df_modifications[
                                                        df_modifications["ModType"]
                                                        == "isoasp"
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                            rx_mapnum = (
                                                df_sequence_overview.loc[
                                                    aa_choice_num, "Mapnum"
                                                ]
                                                + df_sequence_overview.loc[
                                                    aa_choice_num + 1, "Mapnum"
                                                ]
                                            )
                                    else:
                                        if calc_type != "fast":
                                            aa_choice = userinput2_aamod(
                                                aa_options, aamod_finalchoice, path_to_teriyaki
                                            )
                                            if aa_choice == "":
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                            aa_choice_num = int(
                                                [
                                                    part
                                                    for part in aa_choice.split(" ")
                                                    if part
                                                ][0]
                                            )
                                        if (
                                            aamod_finalchoice == "Betaalanine Formation"
                                            or aamod_finalchoice == "bala"
                                        ):
                                            aamod_finalchoice = "Betaalanine Formation"
                                            rx_toperform = (
                                                rdChemReactions.ReactionFromSmarts(
                                                    df_modifications[
                                                        df_modifications["ModType"]
                                                        == "bala"
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                            rx_mapnum = df_sequence_overview.loc[
                                                aa_choice_num, "Mapnum"
                                            ]
                                        elif (
                                            aamod_finalchoice == "Phosphorylation"
                                            or aamod_finalchoice == "phos"
                                        ):
                                            aamod_finalchoice = "Phosphorylation"
                                            rx_toperform = (
                                                rdChemReactions.ReactionFromSmarts(
                                                    df_modifications[
                                                        df_modifications["ModType"]
                                                        == "phos"
                                                    ].iloc[0]["SmartsReaction"]
                                                )
                                            )
                                            rx_mapnum = df_sequence_overview.loc[
                                                aa_choice_num, "Mapnum"
                                            ]

                                rx_mapnum_atomidx = [
                                    atom.GetIdx()
                                    for atom in chain.GetAtoms()
                                    if atom.GetIntProp("mapnum") in rx_mapnum
                                ]

                                for match in chain.GetSubstructMatches(
                                    rx_toperform.GetReactantTemplate(0)
                                ):
                                    if all(idx in rx_mapnum_atomidx for idx in match):
                                        rx_possible = True

                                if rx_possible:

                                    chain_mod, mapnum_store, rx_mapnum = perform_mod_rx(
                                        rx_toperform,
                                        chain_reactant,
                                        additional_reactant,
                                        highest_mapnum_mod,
                                        mapnum_store,
                                        rx_mapnum,
                                        aamod_finalchoice,
                                        df_sequence_overview,
                                        aa_choice_num,
                                        famod_mapnum,
                                    )

                                    if num_equiv_solutions(chain_mod) == 1:
                                        n = 0
                                    else:
                                        if calc_type == "fast":
                                            errormessage = (
                                                "Reaction "
                                                + str(aamod_count + 1)
                                                + " results in "
                                                + str(num_equiv_solutions(chain_mod))
                                                + " products"
                                            )
                                        else:
                                            solution_ignore = tk.messagebox.askyesno(
                                                "Question",
                                                "The reaction would result in "
                                                + str(num_equiv_solutions(chain_mod))
                                                + " products. Do you want to proceed with one of the products?",
                                            )
                                            if solution_ignore:
                                                n = random.randrange(
                                                    0, num_equiv_solutions(chain_mod), 1
                                                )
                                            else:
                                                back = True
                                                aamod_ht_todo -= 1
                                                break
                                    if errormessage == "none":
                                        chain_ext = chain_mod[n][0]
                                        Chem.SanitizeMol(chain_ext)
                                        if (
                                            aamod_finalchoice == "Phosphorylation"
                                            or aamod_finalchoice == "Sulfur Oxidation"
                                        ):
                                            for a in chain_ext.GetAtoms():
                                                if a.HasProp("mapnum") == False:
                                                    highest_mapnum += 1
                                                    a.SetIntProp(
                                                        "mapnum", highest_mapnum
                                                    )
                                                    rx_mapnum.append(highest_mapnum)
                                            mapnum_store.append(highest_mapnum)
                                        chain = chain_ext
                                        df_modif_overview.loc[
                                            aamod_count,
                                            [
                                                "AAModif",
                                                "NumAAModif",
                                                "ModifMyMap",
                                                "ModifType",
                                                "Additional",
                                            ],
                                        ] = [
                                            aa_choice,
                                            aa_choice_num,
                                            rx_mapnum,
                                            aamod_finalchoice,
                                            0,
                                        ]
                                        if (
                                            aamod_finalchoice
                                            == "Fatty Acid Derivatization"
                                        ):
                                            df_modif_overview.loc[
                                                aamod_count, ["Additional"]
                                            ] = fa_info
                                            df_modif_overview.loc[
                                                aamod_count, ["ModifType"]
                                            ] = "Fatty Acid"
                                        elif (
                                            aamod_finalchoice
                                            == "Aspartimide / Isoaspartate"
                                            or aamod_finalchoice == "asp"
                                            or aamod_finalchoice == "isoasp"
                                        ):
                                            df_modif_overview.loc[
                                                aamod_count, ["ModifType"]
                                            ] = option_choice
                                        elif aamod_finalchoice == "Sulfur Oxidation":
                                            df_modif_overview.loc[
                                                aamod_count, ["Additional"]
                                            ] = option_choice[:-6]
                                        elif aamod_finalchoice == "Glycosylation":
                                            df_modif_overview.loc[
                                                aamod_count, ["Additional"]
                                            ] = option_choice
                                        df_sequence_overview.loc[
                                            aa_choice_num, "AminoAcidInfo"
                                        ] += (
                                            " * mod" + "(" + str(aamod_count + 1) + ")"
                                        )
                                        df_sequence_overview.at[
                                            aa_choice_num, "Mapnum"
                                        ] = list(
                                            set(
                                                df_sequence_overview.at[
                                                    aa_choice_num, "Mapnum"
                                                ]
                                                + rx_mapnum
                                            )
                                        )
                                        aamod_count += 1
                                        MyModPeptide = Peptide(
                                            chain,
                                            mapnum_store,
                                            df_sequence_overview,
                                            df_cycl_overview,
                                            df_modif_overview,
                                            assign_pepinfo(chain, df_sequence_overview),
                                        )
                                        peptide_memory.append(MyModPeptide)
                                    else:
                                        window_status, calc_table = errorupdate(
                                            window_status,
                                            calc_type,
                                            calc_table,
                                            fastcalc_pepnum,
                                            errormessage,
                                            new_window_status=2,
                                        )
                                        back = True
                                        aamod = False
                                        if fastcalc_pepnum == len(calc_table):
                                            window_status = 3
                                        else:
                                            break
                                else:
                                    errormessage = aamod_finalchoice + " not possible"
                                    window_status, calc_table = errorupdate(
                                        window_status,
                                        calc_type,
                                        calc_table,
                                        fastcalc_pepnum,
                                        errormessage,
                                        new_window_status=2,
                                    )
                                    aamod = False
                                    back = True
                                    if fastcalc_pepnum == len(calc_table):
                                        window_status = 3
                                    else:
                                        break

                        if calc_type == "fast":
                            if errormessage == "none":
                                df_htpepsummary = create_peptide_summary(peptide_memory)
                                htpepsummary_values = df_htpepsummary.iloc[
                                    :-2, 1
                                ].tolist()
                                calc_table.iloc[
                                    fastcalc_pepnum - 1,
                                    7 : 7 + len(htpepsummary_values),
                                ] = htpepsummary_values
                                calc_table.at[fastcalc_pepnum - 1, "Status"] = "Done"
                            if fastcalc_pepnum == len(calc_table):
                                window_status = 3
                            else:
                                break

                    if window_status == 3:
                        if calc_type == "fast":
                            if (
                                calc_table.iloc[-1, 0]
                                != "Calculation with PICKAPEP (Erckes V. et al. 2024, J. Pept. Sci.)"
                            ):
                                calc_table.loc[fastcalc_pepnum] = [None] * (
                                    len(calc_table.columns)
                                )
                                calc_table.loc[fastcalc_pepnum + 1] = [
                                    "Calculation with PICKAPEP (Erckes V. et al. 2024, J. Pept. Sci.)"
                                ] + [None] * (len(calc_table.columns) - 1)
                                calc_table.loc[
                                    fastcalc_pepnum : fastcalc_pepnum + 1, "Status"
                                ] = ""
                            window_status = userinput3_fastcalc(
                                calc_table, calc_table_path, skip
                            )
                            if window_status == 4:
                                break
                        else:
                            peptide_memory, back, window_status = userinput3(
                                peptide_memory, window_status, path_to_teriyaki, path_to_teriyaki_easteregg
                            )

    except:
        messagebox.showerror(
            "Error", "Sorry, an error occurred. Please restart PICKAPEP."
        )

## ADD DEFAULT PEPTIDE TO DEFAULT PEPTIDE PICKLE FILE ##
# Needed for initialization of default peptide pickle file #

# Keep peptide properties when adding peptide to pickle file
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

# with open('dict_default.pkl', 'rb') as file:
#     dict_default = pickle.load(file)

# print(dict_default)

# dict_default['Somatostatin'] = peptide_memory


# print(dict_default)
# with open('dict_default.pkl', 'wb') as file:
#     pickle.dump(dict_default, file)

# # Check pickle file updates

# angiotensinI_mol = Chem.MolFromSmiles('CC[C@H](C)[C@@H](C(=O)N[C@@H](CC1=CN=CN1)C(=O)N2CCC[C@H]2C(=O)N[C@@H](CC3=CC=CC=C3)C(=O)N[C@@H](CC4=CN=CN4)C(=O)N[C@@H](CC(C)C)C(=O)O)NC(=O)[C@H](CC5=CC=C(C=C5)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(=O)O)N')
# angiotensinII_mol = Chem.MolFromSmiles('CC[C@H](C)[C@@H](C(=O)N[C@@H](CC1=CN=CN1)C(=O)N2CCC[C@H]2C(=O)N[C@@H](CC3=CC=CC=C3)C(=O)O)NC(=O)[C@H](CC4=CC=C(C=C4)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(=O)O)N')
# ciclosporinA_mol = Chem.MolFromSmiles('CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C')
# linaclotide_mol = Chem.MolFromSmiles('C[C@H]1C(=O)N[C@H]2CSSC[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@H](NC(=O)CNC(=O)[C@@H](NC2=O)[C@@H](C)O)C(=O)N[C@@H](CC4=CC=C(C=C4)O)C(=O)O)C(=O)N[C@@H](CSSC[C@@H](C(=O)N3)N)C(=O)N[C@H](C(=O)N5CCC[C@H]5C(=O)N1)CC(=O)N)CC6=CC=C(C=C6)O)CCC(=O)O')
# liraglutide_mol = Chem.MolFromSmiles('CCCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)NCCCC[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C)C(=O)N[C@@H](CC2=CNC3=CC=CC=C32)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)N)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC4=CC=C(C=C4)O)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC5=CC=CC=C5)NC(=O)[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC6=CN=CN6)N)C(=O)O')
# semaglutide_mol = Chem.MolFromSmiles('CC[C@H](C)[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)NC(=O)[C@H](CC3=CC=CC=C3)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@H](C(=O)O)NC(=O)CCCCCCCCCCCCCCCCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)N)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC4=CC=C(C=C4)O)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC5=CC=CC=C5)NC(=O)[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)NC(=O)[C@H](CC6=CN=CN6)N')
# oxytocin_mol = Chem.MolFromSmiles('CC[C@H](C)[C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N1)CC2=CC=C(C=C2)O)N)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N)CC(=O)N)CCC(=O)N')
# desmopressin_mol = Chem.MolFromSmiles('N=C(NCCC[C@H](C(NCC(N)=O)=O)NC([C@@H]1CCCN1C([C@@H]2CSSCCC(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N2)=O)CC(N)=O)=O)CCC(N)=O)=O)Cc3ccccc3)=O)Cc4ccc(O)cc4)=O)=O)=O)N')
# somatostatin_mol = Chem.MolFromSmiles('C[C@H]([C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCCCN)CC2=CNC3=CC=CC=C32)CC4=CC=CC=C4)CC5=CC=CC=C5)CC(=O)N)CCCCN)NC(=O)CNC(=O)[C@H](C)N)C(=O)O)CO)[C@@H](C)O)CC6=CC=CC=C6)O')

# x = semaglutid_mol
# x.HasSubstructMatch(dict_default['Semaglutide'][-1].peptide)

## CREATE PICKAPEP ICON
# filename = r'Teriyaki.png'
# img = Image.open(filename)
# img.save('Teriyaki.ico')

## PYINSTALLER
# Code packing with pyinstaller 
# pyinstaller --onefile --noconsole --add-data="Python_AminoAcids.csv;." --add-data="Python_ModReactions.csv;." --add-data="Python_CyclReactions.csv;." --add-data="Python_TerminalMod.csv;." --add-data="dict_default.pkl;." --add-data="Teriyaki.png;." --add-data="Teriyaki_Easteregg.png;." --add-data="Helpbutton_small.png;." --add-data="Infobutton_small.png;." --add-data="HelpSheet_Userinput1.png;." --add-data="HelpSheet_Userinput1_aminoacids.png;." --add-data="HelpSheet_Userinput1_multiplepepcalc.png;." --add-data="HelpSheet_Userinput2.png;." --add-data="HelpSheet_Userinput2_cycl.png;." --add-data="HelpSheet_Userinput2_mod.png;." --add-data="HelpSheet_Userinput3.png;." --add-data="pymol_path.txt;." --add-data="molecule.sdf;." --splash="Teriyaki_Splash.png" --icon="Teriyaki.ico" --hidden-import="PIL._tkinter_finder" --hidden-import="openpyxl.cell._writer" PICKAPEP.py
