import bibtexparser
import re
import os
import tkinter as tk
from tkinter import IntVar
from tkinter import filedialog
from pylatexenc.latex2text import LatexNodes2Text


def select_file():
    """ファイル選択ダイアログを表示し、選択されたファイルのパスを返します。"""
    return filedialog.askopenfilename()


def parse_bibtex(file_path):
    """指定したBibTexファイルをパースし、bib_databaseを返します。"""
    parser = bibtexparser.bparser.BibTexParser(common_strings=True)

    with open(file_path) as bibtex_file:
        return bibtexparser.load(bibtex_file, parser=parser)


def format_bibliography(entry, index):
    """指定したエントリを書誌情報のフォーマットに整形します。"""
    authors = ', '.join(entry.get('author', '').split(' and '))
    title = convert_special_characters(entry.get('title', ''))  # LaTeXマークアップを処理
    journal = f". {LatexNodes2Text().latex_to_text(entry['journal'])}" if 'journal' in entry else ''
    year = f". {LatexNodes2Text().latex_to_text(entry['year'])}" if 'year' in entry else ''
    volume = f"; {LatexNodes2Text().latex_to_text(entry['volume'])}" if 'volume' in entry else ''
    pages = f": {LatexNodes2Text().latex_to_text(entry['pages'])}" if 'pages' in entry else ''

    return f"{index}. {title}. {authors}{journal}{year}{volume}{pages}\n"

def convert_special_characters(text):
    """LaTeXマークアップを適切に処理します。"""
    text = text.replace(r'$\upalpha$', 'α')  # 特定のマークアップを適切な文字に置換
    text = text.replace(r'$\upbeta$', 'β')  # 特定のマークアップを適切な文字に置換
    text = text.replace(r'$\upgamma$', 'γ')  # 特定のマークアップを適切な文字に置換
    text = text.replace(r'$\updelta$', 'δ')  # 特定のマークアップを適切な文字に置換
    text = text.replace(r'$\upepsilon$', 'ε')  # 特定のマークアップを適切な文字に置換
    text = text.replace(r'$\upzeta$', 'ζ')
    text = text.replace(r'$\uptheta$', 'θ')
    text = text.replace(r'$\upkappa$', 'κ')
    text = text.replace(r'$\uplambda$', 'λ')
    text = text.replace(r'$\upmu$', 'μ')
    text = text.replace(r'$\uppi$', 'π')
    text = text.replace(r'$\upsigma$', 'σ')
    text = text.replace(r'$\uptau$', 'τ')
    text = text.replace(r'$\upphi$', 'φ')
    text = text.replace(r'$\uppsi$', 'ψ')
    text = text.replace(r'$\upomega$', 'ω')
    text = text.replace(r'$\less$sup$\greater$$\mathplus$$\less$/sup$\greater$', '⁺')  # 特定のマークアップを適切な文字に置換
    #text = re.sub(r'\{\\&\}amp\$\\mathsemicolon\$', '&', text)
    #text = text.replace(r'{\&}amp$\mathsemicolon$','&')
    
    text = LatexNodes2Text().latex_to_text(text)  # 一般的なLaTeXマークアップを処理
    text = text.replace(r'&amp', '&')
    return text

def write_bibliography(file_path, bib_database, reverse):
    """書誌情報をフォーマットし、指定したファイルに書き込みます。"""
    directory = os.path.dirname(file_path)
    output_file_path = os.path.join(directory, 'formatted_bibliography.txt')

    with open(output_file_path, 'w', encoding='utf-8') as output_file:
        entries = bib_database.entries
        if reverse.get() == 1:  # チェックがついていれば逆順に
            entries = reversed(entries)
        for i, entry in enumerate(entries, start=1):
            formatted_entry = format_bibliography(entry, i)
            output_file.write(formatted_entry)
    
    print(f"""File was written\n""")


def open_file(reverse_var):
    file_path = select_file()

    if file_path:
        bib_database = parse_bibtex(file_path)
        write_bibliography(file_path, bib_database, reverse_var)


def main():
    root = tk.Tk()
    root.title("BibTex Parser")  # ウィンドウのタイトルを設定
    root.geometry("500x500")
    reverse_var = IntVar()  # チェックボタンの状態を管理する変数
    button = tk.Button(root, text="Open BibTex File", command=lambda: open_file(reverse_var))
    button.pack()

    # チェックボタンの追加
    check = tk.Checkbutton(root, text='Reverse Order', variable=reverse_var)
    check.pack()

    root.mainloop()

if __name__ == "__main__":
    main()