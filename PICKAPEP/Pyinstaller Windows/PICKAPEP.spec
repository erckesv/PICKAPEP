# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['PICKAPEP.py'],
    pathex=[],
    binaries=[],
    datas=[('Python_AminoAcids.csv', '.'), ('Python_ModReactions.csv', '.'), ('Python_CyclReactions.csv', '.'), ('Python_TerminalMod.csv', '.'), ('dict_default.pkl', '.'), ('Teriyaki.png', '.'), ('Teriyaki_Easteregg.png', '.'), ('Helpbutton_small.png', '.'), ('Infobutton_small.png', '.'), ('HelpSheet_Userinput1.png', '.'), ('HelpSheet_Userinput1_aminoacids.png', '.'), ('HelpSheet_Userinput1_multiplepepcalc.png', '.'), ('HelpSheet_Userinput2.png', '.'), ('HelpSheet_Userinput2_cycl.png', '.'), ('HelpSheet_Userinput2_mod.png', '.'), ('HelpSheet_Userinput3.png', '.'), ('pymol_path.txt', '.'), ('molecule.sdf', '.')],
    hiddenimports=['PIL._tkinter_finder', 'openpyxl.cell._writer'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
splash = Splash(
    'Teriyaki_Splash.png',
    binaries=a.binaries,
    datas=a.datas,
    text_pos=None,
    text_size=12,
    minify_script=True,
    always_on_top=True,
)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    splash,
    splash.binaries,
    [],
    name='PICKAPEP',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['Teriyaki.ico'],
)
