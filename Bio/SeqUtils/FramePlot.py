from Bio.Seq import reverse_complement, translate
from Bio.Data import CodonTable


def ascii_frames(seq, table=1, translation_check=True):
    """Return dict with strings representing all possible ORFs in 6-frames for given sequence.

    6 frame translation plot allows explore potential ORFs present depending on chosen coding table.

    The marks are centered on the respective codon. Thus:
           start  codon  stop  stop coding
    mark:    >-    ---     |       -|
    codon:  AUG    TTT    TAG      TAG

    If translation table with ambiguous stop is requested the stop is denoted by "!" instead of "|" and the orf
      translation is not interrupted.

    Frames in dict are as follows:
    keys: [1, 2, 3] for [0, 1, 2] offset respectively in + strand.
    keys: [-1, -2. -3] for [0, 1, 2] offset respectively in - strand.

    Note that this is not an ORF prediction tool!
    Code based on "six_frame_translations" and output appearance inspired by frameplot
      (reference: http://www0.nih.go.jp/~jun/cgi-bin/frameplot.pl, https://doi.org/10.1016/0378-1119(84)90116-1).

    :param seq: [str] sequence to make a plot for
    :param table: [int str CodonTable] genetic code table to use (start codons must be different then stop codons)
    :param translation_check: [bool] if translation of each frame should be conducted
        (raises error when ambiguous bases are present)

    >>> from Bio.SeqUtils.FramePlot import ascii_frames
    >>> res = ascii_frames('ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC')
    >>> for i in [-3, -2, -1, 1, 2, 3]:
    ...    print('"' + res[i] + '"')
    "-------<                         =="
    "=-<                               ="
    "==    |                            "
    "                               | =="
    "=                                 ="
    "== >-----------|              >----"

    >>> res = ascii_frames('ACATGTCAGTCGTGTGATCGTGTACGTCG?ATGATC', translation_check=False)
    >>> for i in [-3, -2, -1, 1, 2, 3]:
    ...    print('"' + res[i] + '"')
    "=-<                               =="
    "==    |                            ="
    "-------<                            "
    "                               >----"
    "=                               | =="
    "== >-----------|                   ="

    >>> res = ascii_frames('ACATGTCAGTCGTGTGATCGTGTACGTCG?ATGATC', translation_check=True)
    Traceback (most recent call last):
    ...
    Bio.Data.CodonTable.TranslationError: Codon 'CG?' is invalid
    """

    try:
        t = CodonTable.ambiguous_dna_by_id[int(table)]
    except ValueError:
        t = CodonTable.ambiguous_dna_by_name[table]
    except (AttributeError, TypeError):
        if isinstance(table, CodonTable.CodonTable):
            t = table
        else:
            raise ValueError('Bad table argument')

    start = set(t.start_codons)
    stop = set(t.stop_codons)

    assert len(start.intersection(stop)) == 0

    # deal with dual stop (remove codons present in dual stop from simple stop)
    dual_stop = {c for c in stop if c in t.forward_table}
    for cod in dual_stop:
        stop.remove(cod)

    anti = reverse_complement(seq)
    length = len(seq)
    desc = {}
    frames = {}
    for i in range(0, 3):
        fragment_length = 3 * ((length - i) // 3)

        desc[i + 1] = _find_orfs(
            seq[i:i + fragment_length],
            start, stop,
            dual_stop
        )
        desc[-(i + 1)] = _find_orfs(
            anti[i:i + fragment_length],
            start, stop,
            dual_stop
        )[::-1].replace('>', '<')

        # check if translation is valid
        if translation_check:
            frames[i + 1] = translate(
                sequence=seq[i:i + fragment_length],
                table=t)
            frames[-(i + 1)] = translate(
                sequence=anti[i:i + fragment_length],
                table=t)[::-1]

    preformat = {}
    lsu = len(seq)
    preformat[3] = '==' + desc[3] + '{}'.format('='*(lsu - len(desc[3]) - 2))
    preformat[2] = '=' + desc[2] + '{}'.format('='*(lsu - len(desc[2]) - 1))
    preformat[1] = desc[1] + '{}'.format('='*(lsu - len(desc[1])))

    p1 = '='*(lsu - len(desc[-1]))
    preformat[-1] = p1 + desc[-1] + '{}'.format('='*(lsu - len(desc[-1]) - len(p1)))
    p2 = '='*(lsu - len(desc[-2]) - 1)
    preformat[-2] = p2 + desc[-2] + '{}'.format('='*(lsu - len(desc[-2]) - len(p2)))
    p3 = '='*(lsu - len(desc[-3]) - 2)
    preformat[-3] = p3 + desc[-3] + '{}'.format('='*(lsu - len(desc[-3]) - len(p3)))

    return preformat


def _find_orfs(seq, start, stop, am_stop):
    """Find and annotate ORF in sequence. Return string representation. The annotation is codon centered.

    >>> ss = 'ACATGTCAGTCGATGTGATCGTGTACGTCGATGATC'
    >>> print(_find_orfs(ss, {'ATG'}, {'TCA'}, {}))
                 >----------------->----

    >>> print(_find_orfs(ss, {'ATG'}, {'TCG'}, {}))
              |  >-----|        |  >----

    >>> print(_find_orfs(ss, {'ATG'}, {}, {'TCG'}))
              !  >-----!--------!-->----

    :param seq: [str] sequence to translate. Len should be multiple of 3.
    :param start: [set list tuple] list of codons to use as start.
    :param stop: [set list tuple] list of codons to use as stop.
    :param am_stop: [set list tuple] list of codons to use as ambiguous stop.

    Note that start stop and am_stop must not include same codons. This is not checked.
    """

    frame = []
    for j in range(0, len(seq), 3):
        cod = seq[j:j + 3]
        if cod in stop:
            frame.append(-1)
        elif cod in start:
            frame.append(1)
        elif cod in am_stop:
            frame.append(-2)
        else:
            frame.append(0)

    c = 0
    orfinit = False
    s = []
    for k in frame:
        if k == 0:
            # normal translation
            if orfinit:
                c += 1
                s.append('---')
            else:
                s.append('   ')
        elif k == 1:
            # start codon
            c += 1
            if orfinit:
                s.append('->-')
            else:
                orfinit = True
                s.append(' >-')
        elif k == -1:
            # stop codon
            c = 0
            if orfinit:
                s.append('-| ')
            else:
                s.append(' | ')
            orfinit = False
        elif k == -2:
            # ambiguous stop codon
            c += 1
            if orfinit:
                s.append('-!-')
            else:
                s.append(' ! ')
        else:
            raise ValueError

    return ''.join(s)


def ascii_frameplot(seq, genetic_code=1, translation_check=True, linelength=60):
    """Return strings representing all possible ORFs in 6-frames for given sequence.

    6 frame translation plot allows explore potential ORFs present depending on chosen coding table.

    The marks are centered on the respective codon. Thus:
           start  codon  stop  stop coding
    mark:    >-    ---     |       -|
    codon:  AUG    TTT    TAG      TAG

    If translation table with ambiguous stop is requested the stop is denoted by "!" instead of "|" and the orf
      translation is not interrupted.

    Note that this is not an ORF prediction tool!
    Code based on "six_frame_translations" and output appearance inspired by frameplot
      (reference: http://www0.nih.go.jp/~jun/cgi-bin/frameplot.pl, https://doi.org/10.1016/0378-1119(84)90116-1).

    :param seq:               [str] sequence to make a plot for
    :param genetic_code:      [int optional] ncbi genetic code to use (start codons must be different then stop codons)
    :param translation_check: [bool optional] if translation of each frame should be conducted
        (raises error when ambiguous bases are present)
    :param linelength:        [int optional] sequence length which is to be printed on one line
        (the actual string is 2 chars longer)

    >>> from Bio.SeqUtils.FramePlot import ascii_frameplot
    >>> res = ascii_frameplot('ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC')
    >>> print(''.join(['"' + i + '"' + chr(10) for i in res.split(chr(10))]))
    "3 == >-----------|              >----"
    "2 =                                 ="
    "1                                | =="
    "5'ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC"
    "3'TGTACAGTCAGCACACTAGCACATGCAGCTACTAG"
    "1 ==    |                            "
    "2 =-<                               ="
    "3 -------<                         =="
    ""
    ""
    <BLANKLINE>

    Use different translation table.

    >>> res = ascii_frameplot('ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC', 11)
    >>> print(''.join(['"' + i + '"' + chr(10) for i in res.split(chr(10))]))
    "3 == >-------->--|              >-->-"
    "2 =             >-->-->-------------="
    "1                                | =="
    "5'ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC"
    "3'TGTACAGTCAGCACACTAGCACATGCAGCTACTAG"
    "1 ==    |                            "
    "2 =-<--------------------------<--< ="
    "3 -------<--------<                =="
    ""
    ""
    <BLANKLINE>

    Use ambiguous translation table. Note that BiopythonWarning is also raised.

    >>> res = ascii_frameplot('ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC', 27)
    >>> print(''.join(['"' + i + '"' + chr(10) for i in res.split(chr(10))]))
    "3 == >-----------!-------------->----"
    "2 =                                 ="
    "1                                ! =="
    "5'ACATGTCAGTCGTGTGATCGTGTACGTCGATGATC"
    "3'TGTACAGTCAGCACACTAGCACATGCAGCTACTAG"
    "1 ==    !                            "
    "2 =-<                               ="
    "3                                  =="
    ""
    ""
    <BLANKLINE>

    Show behavior with non-standard bases.

    >>> res = ascii_frameplot('ACATGTCAGTCGTGTGATCGTTGAGTACGTCGAT?ATC', translation_check=False)
    >>> print(''.join(['"' + i + '"' + chr(10) for i in res.split(chr(10))]))
    "3 == >-----------|     >----------------"
    "2 =                                    ="
    "1                       |             =="
    "5'ACATGTCAGTCGTGTGATCGTTGAGTACGTCGAT?ATC"
    "3'TGTACAGTCAGCACACTAGCAACTCATGCAGCTA?TAG"
    "1 ==    |                               "
    "2 =-<                                  ="
    "3 -------<                            =="
    ""
    ""
    <BLANKLINE>

    >>> ascii_frameplot('ACATGTCAGTCGTGTGATCGTTGAGTACGTCGAT?ATC', translation_check=True)
    Traceback (most recent call last):
    ...
    Bio.Data.CodonTable.TranslationError: Codon 'T?A' is invalid


    >>> res = ascii_frameplot('CGATCGTAGTCGATCGTGACTAGTCAGATCGCTAGTCGTAGCTACGTACGTCGATGC', 11, linelength=30)
    >>> print(''.join(['"' + i + '"' + chr(10) for i in res.split(chr(10))]))
    "3 == >-----------------|        "
    "2 =                |            "
    "1        |     >-->----------->-"
    "5'CGATCGTAGTCGATCGTGACTAGTCAGATC"
    "3'GCTAGCATCAGCTAGCACTGATCAGTCTAG"
    "1 -------------------------<    "
    "2 ==----------<           |--<  "
    "3 =-<                 |         "
    ""
    "3    |     |              >-="
    "2                          =="
    "1 ---------------------------"
    "5'GCTAGTCGTAGCTACGTACGTCGATGC"
    "3'CGATCAGCATCGATGCATGCAGCTACG"
    "1                            "
    "2             |             ="
    "3   |--------------------< =="
    ""
    ""
    <BLANKLINE>
    """
    desc = ascii_frames(seq, genetic_code, translation_check)

    comp = reverse_complement(seq)[::-1]

    res = ''
    for i in range(0, len(seq), linelength):
        subseq = seq[i:i + linelength]
        csubseq = comp[i:i + linelength]

        res += '3 ' + desc[3][i:i + linelength] + '\n'
        res += '2 ' + desc[2][i:i + linelength] + '\n'
        res += '1 ' + desc[1][i:i + linelength] + '\n'

        res += "5'" + subseq + '\n'
        res += "3'" + csubseq + '\n'

        res += '1 ' + desc[-1][i:i + linelength] + '\n'
        res += '2 ' + desc[-2][i:i + linelength] + '\n'
        res += '3 ' + desc[-3][i:i + linelength] + '\n\n'
    return res

if __name__ == "__main__":
    import doctest
    from Bio._utils import run_doctest
    run_doctest(optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
