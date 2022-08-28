#!/bin/python3

from pathlib import Path
import pandas as pd
import csv, os

from tkinter import *
from tkinter import ttk, filedialog

from MaxQuantHandler.mq_utils import runner_utils as ru
from filter_ids import filter_protein_ids as fpi
from get_uniprot_mapping import get_uniprot_mappings as gum
from remap_genenames import remap_genenames as rg


OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path("./assets")


def maxquant_handler():
    home()


def home():
    window = Tk()
    window.title("MaxQuant Handler - Home")
    window.geometry("1043x546")
    window.configure(bg="#008E85")

    canvas = Canvas(
        window,
        bg="#008E85",
        height=546,
        width=1043,
        bd=0,
        highlightthickness=0,
        relief="ridge"
    )

    canvas.place(x=0, y=0)
    canvas.create_text(
        316.0,
        20.0,
        anchor="nw",
        text="MaxQuant Handler",
        fill="#FFFFFF",
        font=("Georgia", 48 * -1)
    )

    canvas.create_rectangle(
        44.0,
        135.0,
        331.0,
        511.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_rectangle(
        706.0,
        135.0,
        993.0,
        511.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_rectangle(
        375.0,
        135.0,
        662.0,
        511.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_text(
        245.0,
        82.0,
        anchor="nw",
        text="Python-based helpers to clean, filter and fill your MaxQuant file.",
        fill="#FFFFFF",
        font=("Georgia", 20 * -1)
    )

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: filter_protein_ids(tk_window=window),
        relief="flat"
    )
    button_1.place(
        x=107.0,
        y=449.0,
        width=160.0,
        height=43.0
    )

    button_image_2 = PhotoImage(
        file=relative_to_assets("button_2.png"))
    button_2 = Button(
        image=button_image_2,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: get_uniprot_mapping(tk_window=window),
        relief="flat"
    )
    button_2.place(
        x=770.0,
        y=449.0,
        width=160.0,
        height=43.0
    )

    button_image_3 = PhotoImage(
        file=relative_to_assets("button_3.png"))
    button_3 = Button(
        image=button_image_3,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: remap_genenames(tk_window=window),
        relief="flat"
    )
    button_3.place(
        x=439.0,
        y=449.0,
        width=160.0,
        height=43.0
    )

    canvas.create_text(
        96.0,
        165.0,
        anchor="nw",
        text="Filter Protein IDs",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_text(
        736.0,
        165.0,
        anchor="nw",
        text="Get Uniprot Mapping",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_text(
        418.0,
        165.0,
        anchor="nw",
        text="Remap Gene Names",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_text(
        73.0,
        269.0,
        anchor="nw",
        text="Filter the Protein IDs in your",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        73.0,
        290.0,
        anchor="nw",
        text="MaxQuant File based on",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        73.0,
        311.0,
        anchor="nw",
        text="Organism. Optionally you",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        73.0,
        332.0,
        anchor="nw",
        text="can also input a file with a",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        73.0,
        353.0,
        anchor="nw",
        text="single column.",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        259.0,
        anchor="nw",
        text="Remap the Gene Names in the",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        280.0,
        anchor="nw",
        text="MaxQuant File using either a",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        301.0,
        anchor="nw",
        text="given Fasta File, Mappings",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        322.0,
        anchor="nw",
        text="from the Uniprot database or",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        343.0,
        anchor="nw",
        text="both. Optionally you can also",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        364.0,
        anchor="nw",
        text="simply get the Header",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        385.0,
        anchor="nw",
        text="Information from the Fasta",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        404.0,
        406.0,
        anchor="nw",
        text="File.",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        259.0,
        anchor="nw",
        text="Get a complete mapping from",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        280.0,
        anchor="nw",
        text="UniProt based on either given",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        301.0,
        anchor="nw",
        text="ProteinIDs or Gene Names of",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        322.0,
        anchor="nw",
        text="the MaxQuant File with",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        343.0,
        anchor="nw",
        text="Information about Organism,",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        364.0,
        anchor="nw",
        text="Review Status associated IDs.",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        385.0,
        anchor="nw",
        text="Optionally you can also input",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        735.0,
        406.0,
        anchor="nw",
        text="a file with a single column.",
        fill="#646464",
        font=("Georgia", 16 * -1)
    )
    window.resizable(False, False)
    window.mainloop()


def filter_protein_ids(tk_window: Tk):
    tk_window.destroy()
    window = Tk()
    window.title("MaxQuant Handler - Filter Protein IDs")
    window.geometry("1043x546")
    window.configure(bg="#008E85")

    canvas = Canvas(
        window,
        bg="#008E85",
        height=546,
        width=1043,
        bd=0,
        highlightthickness=0,
        relief="ridge"
    )

    canvas.place(x=0, y=0)
    canvas.create_text(
        316.0,
        20.0,
        anchor="nw",
        text="MaxQuant Handler",
        fill="#FFFFFF",
        font=("Georgia", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    def run_script():
        if outdir_path.get().startswith("Set output"):
            outdir_path.set("./")
        data = pd.read_table(data_path.get(), sep=ru.find_delimiter(data_path.get())).fillna("")
        file_name = Path(data_path.get()).stem
        df = fpi(data = data, organism = ru.organisms[organism.get()], decoy = decoy.get()==1,
                 action = action.get(), reviewed = reviewed.get()==1)
        df.to_csv(os.path.join(outdir_path.get(), file_name + "_filtered.txt"), header=True, index=False,
                  quoting=csv.QUOTE_NONNUMERIC, sep=" ")

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=run_script,
        relief="flat"
    )
    button_1.place(
        x=703.0,
        y=452.0,
        width=160.0,
        height=43.0
    )
    combostyle = ttk.Style()
    combostyle.theme_create('combostyle', parent='alt', settings={'TCombobox': {'configure':
                                                                                    {'selectbackground': 'white',
                                                                                     'fieldbackground': 'white',
                                                                                     'selectforeground': 'black',
                                                                                     'background': 'white'}}})
    combostyle.theme_use('combostyle')
    canvas.create_text(
        734.0,
        117.0,
        anchor="nw",
        text="Required ",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )
    ttk.Label(window, text="Select your MaxQuant File from your computer.", background="white").place(
        x=600.0,
        y=150.0,
        height=30.0
    )
    data_path = StringVar()
    data_path.set("Path to data file ...")
    def select_file():
        data_path.set(filedialog.askopenfilename())
    Entry(master=window, textvariable=data_path).place(
        x=603.0,
        y=185.0,
        width=280.0,
        height=33.0
    )
    outdir_button = Button(window, text="Browse", command=select_file, background="#FFFFFF")
    outdir_button.place(
        x=880.0,
        y=180.0,
        width=80.0,
        height=43.0
    )
    ttk.Label(window, text="Select organism:", background="white").place(
        x=600.0,
        y=230.0,
        height=30.0
    )
    organism = StringVar()
    organism_box = ttk.Combobox(window, width=27, textvariable=organism)
    organism_box['state'] = 'readonly'
    organism_box['values'] = ('human', 'rat', 'mouse')
    organism_box.current(0)
    organism_box.place(
        x=810.0,
        y=230.0,
        width=150.0,
        height=30.0
    )
    ttk.Label(window, text="Select action for empty cells:", background="white").place(
        x=600.0,
        y=310.0,
        height=30.0
    )
    action = StringVar()
    action_box = ttk.Combobox(window, width=27, textvariable=action)
    action_box['state'] = 'readonly'
    action_box['values'] = ('delete', 'fill', 'keep')
    action_box.current(0)
    action_box.place(
        x=810.0,
        y=310.0,
        width=150.0,
        height=30.0
    )
    canvas.create_rectangle(
        581.0,
        131.0,
        712.0,
        133.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        131.0,
        984.0,
        133.0,
        fill="#646464",
        outline="")

    canvas.create_text(
        734.0,
        277.0,
        anchor="nw",
        text="Optional",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    decoy = IntVar()
    Checkbutton(window, text="Keep decoy IDs", variable=decoy, bg="white").place(
        x=810.0,
        y=350.0,
        width=150.0,
        height=43.0
    )

    reviewed = IntVar()
    Checkbutton(window, text="Take reviewed IDs", variable=reviewed, bg="white").place(
        x=600.0,
        y=350.0,
        width=150.0,
        height=43.0
    )

    outdir_path = StringVar()
    outdir_path.set("Set output directory. Default: ./")
    def select_dir():
        outdir_path.set( filedialog.askdirectory())
    Entry(master=window, textvariable=outdir_path).place(
        x=603.0,
        y=405.0,
        width=280.0,
        height=33.0
    )
    outdir_button = Button(window, text="Browse", command=select_dir, background="#FFFFFF")
    outdir_button.place(
        x=880.0,
        y=400.0,
        width=80.0,
        height=43.0
    )


    canvas.create_text(
        148.0,
        118.0,
        anchor="nw",
        text="Filter Protein IDs",
        fill="#FFFFFF",
        font=("Georgia", 22 * -1)
    )

    canvas.create_rectangle(
        581.0,
        292.0,
        712.0,
        294.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        292.0,
        984.0,
        294.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        0.0,
        104.0,
        552.0,
        109.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_text(
        33.0,
        208.0,
        anchor="nw",
        text="“Protein IDs”",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        240.0,
        anchor="nw",
        text="What to do, if a cell is empty after filtering:",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        261.0,
        anchor="nw",
        text="keep: keep the empty row",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        282.0,
        anchor="nw",
        text="fill: try to map Protein IDs from Gene Name",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        303.0,
        anchor="nw",
        text="if column with “Gene names” is given",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        324.0,
        anchor="nw",
        text="delete: delete the row",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        166.0,
        anchor="nw",
        text="Select the File with the Protein IDs. Please make sure that the",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="column is called:",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        365.0,
        anchor="nw",
        text="Set checkmark on reviewed, if during the fill option only",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        386.0,
        anchor="nw",
        text="reviewed Protein IDs should be chosen.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        438.0,
        anchor="nw",
        text="Set checkmark on decoy, if IDs with REV_, CON_ should be",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        459.0,
        anchor="nw",
        text="deleted.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )
    window.resizable(False, False)
    window.mainloop()


def remap_genenames(tk_window: Tk):
    tk_window.destroy()
    window = Tk()
    window.title("MaxQuant Handler - Remap Gene Names")
    window.geometry("1043x546")
    window.configure(bg="#008E85")

    canvas = Canvas(
        window,
        bg="#008E85",
        height=546,
        width=1043,
        bd=0,
        highlightthickness=0,
        relief="ridge"
    )

    canvas.place(x=0, y=0)
    canvas.create_text(
        316.0,
        20.0,
        anchor="nw",
        text="MaxQuant Handler",
        fill="#FFFFFF",
        font=("Georgia", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    def run_script():
        if outdir_path.get().startswith("Set output"):
            outdir_path.set("./")
        data = pd.read_table(data_path.get(), sep=ru.find_delimiter(data_path.get())).fillna("")
        file_name = Path(data_path.get()).stem
        df = rg(data=data, mode=mode.get(), skip_filled=fill.get()==1, organism=ru.organisms[organism.get()],
                fasta=fasta.get())
        df.to_csv(os.path.join(outdir_path.get(), file_name + "_filtered.txt"), header=True, index=False,
                  quoting=csv.QUOTE_NONNUMERIC, sep=" ")

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=run_script,
        relief="flat"
    )
    button_1.place(
        x=703.0,
        y=452.0,
        width=160.0,
        height=43.0
    )

    canvas.create_text(
        734.0,
        117.0,
        anchor="nw",
        text="Required ",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_rectangle(
        581.0,
        131.0,
        712.0,
        133.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        131.0,
        984.0,
        133.0,
        fill="#646464",
        outline="")

    ttk.Label(window, text="Select your MaxQuant File from your computer.", background="white").place(
        x=600.0,
        y=150.0,
        height=30.0
    )
    data_path = StringVar()
    data_path.set("Path to data file ...")
    def select_file():
        data_path.set(filedialog.askopenfilename())

    combostyle = ttk.Style()
    combostyle.theme_create('combostyle', parent='alt', settings={'TCombobox': {'configure':
                                                                                    {'selectbackground': 'white',
                                                                                     'fieldbackground': 'white',
                                                                                     'selectforeground': 'black',
                                                                                     'background': 'white'}}})
    combostyle.theme_use('combostyle')
    ttk.Label(window, text="Select mode (see description):", background="white").place(
        x=600.0,
        y=230.0,
        height=30.0
    )
    mode = StringVar()
    mode_box = ttk.Combobox(window, width=27, textvariable=mode)
    mode_box['state'] = 'readonly'
    mode_box['values'] = ('all', 'fasta', 'uniprot', 'uniprot_one')
    mode_box.current(0)
    mode_box.place(
        x=810.0,
        y=230.0,
        width=150.0,
        height=30.0
    )

    Entry(master=window, textvariable=data_path).place(
        x=603.0,
        y=185.0,
        width=280.0,
        height=33.0
    )
    file_button = Button(window, text="Browse", command=select_file, background="#FFFFFF")
    file_button.place(
        x=880.0,
        y=180.0,
        width=80.0,
        height=43.0
    )

    canvas.create_text(
        734.0,
        277.0,
        anchor="nw",
        text="Optional ",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_text(
        174.0,
        123.0,
        anchor="nw",
        text="Remap Genenames",
        fill="#FFFFFF",
        font=("Georgia", 22 * -1)
    )

    canvas.create_rectangle(
        581.0,
        292.0,
        712.0,
        294.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        292.0,
        984.0,
        294.0,
        fill="#646464",
        outline="")

    fill = IntVar()
    Checkbutton(window, text="Fill", variable=fill, bg="white").place(
        x=600.0,
        y=310.0,
        width=60.0,
        height=43.0
    )

    fasta = StringVar()
    fasta.set("Optionally select FASTA file.")

    def select_fasta():
        fasta.set(filedialog.askopenfilename())

    Entry(master=window, textvariable=fasta).place(
        x=683.0,
        y=315.0,
        width=200.0,
        height=33.0
    )
    fasta_button = Button(window, text="Browse", command=select_fasta, background="#FFFFFF")
    fasta_button.place(
        x=880.0,
        y=310.0,
        width=80.0,
        height=43.0
    )

    ttk.Label(window, text="Select organism:", background="white").place(
        x=600.0,
        y=360.0,
        height=30.0
    )
    organism = StringVar()
    organism_box = ttk.Combobox(window, width=27, textvariable=organism)
    organism_box['state'] = 'readonly'
    organism_box['values'] = ('human', 'rat', 'mouse')
    organism_box.current(0)
    organism_box.place(
        x=810.0,
        y=360.0,
        width=150.0,
        height=30.0
    )

    outdir_path = StringVar()
    outdir_path.set("Set output directory. Default: ./")

    def select_dir():
        outdir_path.set(filedialog.askdirectory())

    Entry(master=window, textvariable=outdir_path).place(
        x=603.0,
        y=405.0,
        width=280.0,
        height=33.0
    )
    outdir_button = Button(window, text="Browse", command=select_dir, background="#FFFFFF")
    outdir_button.place(
        x=880.0,
        y=400.0,
        width=80.0,
        height=43.0
    )
    canvas.create_rectangle(
        0.0,
        104.0,
        552.0,
        109.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_text(
        33.0,
        208.0,
        anchor="nw",
        text="Fasta File if it should be used in the chosen Modes.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        250.0,
        anchor="nw",
        text="all: use primarly fasta infos and additionally uniprot infos.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        271.0,
        anchor="nw",
        text="fasta: use information extracted from fasta headers.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        292.0,
        anchor="nw",
        text="uniprot:	use mapping information from uniprot and use all",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        313.0,
        anchor="nw",
        text="gene names.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        402.0,
        anchor="nw",
        text="Set checkmark on fill, if only rows without previously mapped",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        423.0,
        anchor="nw",
        text="Gene Names should be mapped.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        460.0,
        anchor="nw",
        text="Optionally select an Organism, so only Organism associated",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        481.0,
        anchor="nw",
        text="Protein IDs will be used during the Mapping.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        334.0,
        anchor="nw",
        text="uniprot_one: use mapping information from uniprot and only",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        355.0,
        anchor="nw",
        text="use most frequent single gene name.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        166.0,
        anchor="nw",
        text="Remap Gene Names in a MaxQuant File based on one of four",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="possible Modes. Keep in mind, that you have to provide the",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    window.resizable(False, False)
    window.mainloop()


def get_uniprot_mapping(tk_window: Tk):
    tk_window.destroy()
    window = Tk()
    window.title("MaxQuant Handler - Get UniProt Mapping")
    window.geometry("1043x546")
    window.configure(bg="#008E85")

    canvas = Canvas(
        window,
        bg="#008E85",
        height=546,
        width=1043,
        bd=0,
        highlightthickness=0,
        relief="ridge"
    )

    canvas.place(x=0, y=0)
    canvas.create_text(
        316.0,
        20.0,
        anchor="nw",
        text="MaxQuant Handler",
        fill="#FFFFFF",
        font=("Georgia", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    def run_script():
        if outdir_path.get().startswith("Set output"):
            outdir_path.set("./")
        data = pd.read_table(data_path.get(), sep=ru.find_delimiter(data_path.get())).fillna("")
        file_name = Path(data_path.get()).stem
        df = gum(data=data, in_type = input_type.get(),organism=ru.organisms[organism.get()])
        df.to_csv(os.path.join(outdir_path.get(), file_name + "_uniprot_" + input_type.get() + "_mapping.csv"),
                  header=True, index=False)

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=run_script,
        relief="flat"
    )
    button_1.place(
        x=703.0,
        y=452.0,
        width=160.0,
        height=43.0
    )

    canvas.create_text(
        734.0,
        117.0,
        anchor="nw",
        text="Required ",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )

    canvas.create_rectangle(
        581.0,
        131.0,
        712.0,
        133.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        131.0,
        984.0,
        133.0,
        fill="#646464",
        outline="")
    ttk.Label(window, text="Select your MaxQuant File from your computer.", background="white").place(
        x=600.0,
        y=150.0,
        height=30.0
    )
    data_path = StringVar()
    data_path.set("Path to data file ...")

    def select_file():
        data_path.set(filedialog.askopenfilename())

    Entry(master=window, textvariable=data_path).place(
        x=603.0,
        y=185.0,
        width=280.0,
        height=33.0
    )
    file_button = Button(window, text="Browse", command=select_file, background="#FFFFFF")
    file_button.place(
        x=880.0,
        y=180.0,
        width=80.0,
        height=43.0
    )

    combostyle = ttk.Style()
    combostyle.theme_create('combostyle', parent='alt', settings={'TCombobox': {'configure':
                                                                                    {'selectbackground': 'white',
                                                                                     'fieldbackground': 'white',
                                                                                     'selectforeground': 'black',
                                                                                     'background': 'white'}}})
    combostyle.theme_use('combostyle')
    ttk.Label(window, text="Select organism:", background="white").place(
        x=600.0,
        y=230.0,
        height=30.0
    )
    organism = StringVar()
    organism_box = ttk.Combobox(window, width=27, textvariable=organism)
    organism_box['state'] = 'readonly'
    organism_box['values'] = ('human', 'rat', 'mouse')
    organism_box.current(0)
    organism_box.place(
        x=810.0,
        y=230.0,
        width=150.0,
        height=30.0
    )
    ttk.Label(window, text="Select you Input type:", background="white").place(
        x=600.0,
        y=280.0,
        height=30.0
    )
    input_type = StringVar()
    input_type_box = ttk.Combobox(window, width=27, textvariable=input_type)
    input_type_box['state'] = 'readonly'
    input_type_box['values'] = ('proteinID', 'genename')
    input_type_box.current(0)
    input_type_box.place(
        x=810.0,
        y=280.0,
        width=150.0,
        height=30.0
    )
    canvas.create_text(
        734.0,
        327.0,
        anchor="nw",
        text="Optional ",
        fill="#646464",
        font=("Georgia", 22 * -1)
    )
    canvas.create_rectangle(
        581.0,
        342.0,
        712.0,
        344.0,
        fill="#646464",
        outline="")

    canvas.create_rectangle(
        853.0,
        342.0,
        984.0,
        344.0,
        fill="#646464",
        outline="")
    ttk.Label(window, text="Select your output directory.", background="white").place(
        x=600.0,
        y=370.0,
        height=30.0
    )
    outdir_path = StringVar()
    outdir_path.set("Set output directory. Default: ./")

    def select_dir():
        outdir_path.set(filedialog.askdirectory())

    Entry(master=window, textvariable=outdir_path).place(
        x=603.0,
        y=405.0,
        width=280.0,
        height=33.0
    )
    outdir_button = Button(window, text="Browse", command=select_dir, background="#FFFFFF")
    outdir_button.place(
        x=880.0,
        y=400.0,
        width=80.0,
        height=43.0
    )
    canvas.create_text(
        159.0,
        123.0,
        anchor="nw",
        text="Get UniProt Mapping",
        fill="#FFFFFF",
        font=("Georgia", 22 * -1)
    )

    canvas.create_rectangle(
        0.0,
        104.0,
        552.0,
        109.0,
        fill="#FFFFFF",
        outline="")

    canvas.create_text(
        33.0,
        166.0,
        anchor="nw",
        text="Get the UniProt Mapping for either Protein IDs or Gene",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="Names.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        208.0,
        anchor="nw",
        text="MaxQuant Files and Files with a single Column are supported.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        262.0,
        anchor="nw",
        text="Please make sure the columns have the correct name:",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        283.0,
        anchor="nw",
        text="“Protein IDs” for Protein IDs",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        304.0,
        anchor="nw",
        text="“Gene names” for Gene Names",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        358.0,
        anchor="nw",
        text="Select your Input Type:",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        379.0,
        anchor="nw",
        text="“proteinID” for Protein IDs",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        400.0,
        anchor="nw",
        text="“genename” for Gene Names",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        452.0,
        anchor="nw",
        text="Lastly, select the desired Organism. This is important, as the",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )

    canvas.create_text(
        33.0,
        473.0,
        anchor="nw",
        text="Database has too many Organisms to get the IDs from all.",
        fill="#FFFFFF",
        font=("Georgia", 16 * -1)
    )
    window.resizable(False, False)
    window.mainloop()


def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)


if __name__ == "__main__":
    maxquant_handler()
