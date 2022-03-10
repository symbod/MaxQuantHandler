#!/bin/python3

from pathlib import Path

from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage

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
        font=("Poly", 48 * -1)
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
        font=("WorkSans Regular", 18 * -1)
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
        font=("Poly", 24 * -1)
    )

    canvas.create_text(
        736.0,
        165.0,
        anchor="nw",
        text="Get Uniprot Mapping",
        fill="#646464",
        font=("Poly", 24 * -1)
    )

    canvas.create_text(
        418.0,
        165.0,
        anchor="nw",
        text="Remap Gene Names",
        fill="#646464",
        font=("Poly", 24 * -1)
    )

    canvas.create_text(
        73.0,
        269.0,
        anchor="nw",
        text="Filter the Protein IDs in your",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        73.0,
        290.0,
        anchor="nw",
        text="MaxQuant File based on",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        73.0,
        311.0,
        anchor="nw",
        text="Organism. Optionally you",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        73.0,
        332.0,
        anchor="nw",
        text="can also input a file with a",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        73.0,
        353.0,
        anchor="nw",
        text="single column.",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        259.0,
        anchor="nw",
        text="Remap the Gene Names in the",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        280.0,
        anchor="nw",
        text="MaxQuant File using either a",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        301.0,
        anchor="nw",
        text="given Fasta File, Mappings",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        322.0,
        anchor="nw",
        text="from the Uniprot database or",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        343.0,
        anchor="nw",
        text="both. Optionally you can also",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        364.0,
        anchor="nw",
        text="simply get the Header",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        385.0,
        anchor="nw",
        text="Information from the Fasta",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        404.0,
        406.0,
        anchor="nw",
        text="File.",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        259.0,
        anchor="nw",
        text="Get a complete mapping from",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        280.0,
        anchor="nw",
        text="UniProt based on either given",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        301.0,
        anchor="nw",
        text="ProteinIDs or Gene Names of",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        322.0,
        anchor="nw",
        text="the MaxQuant File with",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        343.0,
        anchor="nw",
        text="Information about Organism,",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        364.0,
        anchor="nw",
        text="Review Status associated IDs.",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        385.0,
        anchor="nw",
        text="Optionally you can also input",
        fill="#646464",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        735.0,
        406.0,
        anchor="nw",
        text="a file with a single column.",
        fill="#646464",
        font=("Poly", 18 * -1)
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
        font=("Poly", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: print("button_1 clicked"),
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
        font=("Poly", 24 * -1)
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
        text="Optional ",
        fill="#646464",
        font=("Poly", 24 * -1)
    )

    canvas.create_text(
        148.0,
        118.0,
        anchor="nw",
        text="Filter Protein IDs",
        fill="#FFFFFF",
        font=("Poly", 24 * -1)
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
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        240.0,
        anchor="nw",
        text="What to do, if a cell is empty after filtering:",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        261.0,
        anchor="nw",
        text="keep: keep the empty row",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        282.0,
        anchor="nw",
        text="fill: try to map Protein IDs from Gene Name",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        303.0,
        anchor="nw",
        text="if column with “Gene names” is given",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        324.0,
        anchor="nw",
        text="delete: delete the row",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        166.0,
        anchor="nw",
        text="Select the File with the Protein IDs. Please make sure that the",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="column is called:",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        365.0,
        anchor="nw",
        text="Set checkmark on reviewed, if during the fill option only",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        386.0,
        anchor="nw",
        text="reviewed Protein IDs should be chosen.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        438.0,
        anchor="nw",
        text="Set checkmark on decoy, if IDs with REV_, CON_ should be",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        459.0,
        anchor="nw",
        text="deleted.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
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
        font=("Poly", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: print("button_1 clicked"),
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
        font=("Poly", 24 * -1)
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
        text="Optional ",
        fill="#646464",
        font=("Poly", 24 * -1)
    )

    canvas.create_text(
        174.0,
        123.0,
        anchor="nw",
        text="Remap Genenames",
        fill="#FFFFFF",
        font=("Poly", 24 * -1)
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
        text="Fasta File if it should be used in the chosen Modes.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        250.0,
        anchor="nw",
        text="all: use primarly fasta infos and additionally uniprot infos.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        271.0,
        anchor="nw",
        text="fasta: use information extracted from fasta headers.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        292.0,
        anchor="nw",
        text="uniprot:	use mapping information from uniprot and use all",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        313.0,
        anchor="nw",
        text="gene names.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        402.0,
        anchor="nw",
        text="Set checkmark on fill, if only rows without previously mapped",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        423.0,
        anchor="nw",
        text="Gene Names should be mapped.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        460.0,
        anchor="nw",
        text="Optionally select an Organism, so only Organism associated",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        481.0,
        anchor="nw",
        text="Protein IDs will be used during the Mapping.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        334.0,
        anchor="nw",
        text="uniprot_one: use mapping information from uniprot and only",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        355.0,
        anchor="nw",
        text="use most frequent single gene name.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        166.0,
        anchor="nw",
        text="Remap Gene Names in a MaxQuant File based on one of four",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="possible Modes. Keep in mind, that you have to provide the",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
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
        font=("Poly", 48 * -1)
    )

    canvas.create_rectangle(
        552.0,
        104.0,
        1014.0,
        522.0,
        fill="#FFFFFF",
        outline="")

    button_image_1 = PhotoImage(
        file=relative_to_assets("button_1.png"))
    button_1 = Button(
        image=button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=lambda: print("button_1 clicked"),
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
        font=("Poly", 24 * -1)
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
        159.0,
        123.0,
        anchor="nw",
        text="Get UniProt Mapping",
        fill="#FFFFFF",
        font=("Poly", 24 * -1)
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
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        187.0,
        anchor="nw",
        text="Names.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        208.0,
        anchor="nw",
        text="MaxQuant Files and Files with a single Column are supported.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        262.0,
        anchor="nw",
        text="Please make sure the columns have the correct name:",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        283.0,
        anchor="nw",
        text="“Protein IDs” for Protein IDs",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        304.0,
        anchor="nw",
        text="“Gene names” for Gene Names",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        358.0,
        anchor="nw",
        text="Select your Input Type:",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        379.0,
        anchor="nw",
        text="“proteinID” for Protein IDs",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        400.0,
        anchor="nw",
        text="“genename” for Gene Names",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        452.0,
        anchor="nw",
        text="Lastly, select the desired Organism. This is important, as the",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )

    canvas.create_text(
        33.0,
        473.0,
        anchor="nw",
        text="Database has too many Organisms to get the IDs from all.",
        fill="#FFFFFF",
        font=("Poly", 18 * -1)
    )
    window.resizable(False, False)
    window.mainloop()


def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)


if __name__ == "__main__":
    maxquant_handler()
