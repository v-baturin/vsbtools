import json
import re


def picklines(thefile, whatlines):
    return [x for i, x in enumerate(thefile) if i in whatlines]


def getblock(thefile, start_line, end_line):
    # Gets the LAST text block between lines containing 'start_line' and 'end_line'
    start_no = end_no = start_no_raw = end_no_raw = -1
    with open(thefile) as logfile:
        for num, line in enumerate(logfile):
            if start_line in line:
                start_no_raw = num
            elif end_line in line:
                end_no_raw = num
            if end_no_raw > start_no_raw:
                start_no = start_no_raw
                end_no = end_no_raw
        if start_no == -1 or end_no == -1:
            logfile.seek(0, 0)
            print("Cannot find full block: " + thefile + ": Empty list returned")
            return []
        logfile.seek(0, 0)
        return picklines(logfile, range(start_no + 1, end_no))


def gauformat2dict(string: str, gau_route=False):
    out_dict = dict()
    string = string.casefold().strip()
    string = re.sub(' +', ' ', string)
    string = re.sub('\s*=\s*', '=', string)
    string = re.sub('\s*\)', ')', string)
    string = re.sub('\s*\(\s*', '(', string)
    string = re.sub('\s*,\s*', ',', string)
    string = re.sub(r'([^=])\(', r'\1=(', string)

    if gau_route:
        string = re.sub(r'^#\w? +', '', string)
        list_from_string = string.split()
        out_dict.update({'approach_basis': re.sub('^#', '', list_from_string[0])})
        del list_from_string[0]
    else:
        list_from_string = string.split()
    step = ','.join(list_from_string)
    step = re.sub(r'([^ \(\),=]+)', r"'\1'", step)
    step = re.sub('([,\(])(\'[^ \(\),=]+\')(?=([,\)]|$))', r"\1\2:''",
                  step)  # insertion of colons for params without options
    step = re.sub('=\(', ':{', step)
    step = re.sub('\)', '}', step)
    finalstr = '{' + re.sub('=', ':', step) + '}'
    out_dict.update(json.loads(finalstr.replace('\'', '\"')))
    return out_dict


def recursive_print(dct: dict):
    txt1 = ''
    for k, v in dct.items():
        if isinstance(v, dict):
            txt1 += (k + '=(' + recursive_print(v) + ')')
        elif bool(v):
            txt1 += (k + '=' + v + ',')
        else:
            txt1 += (k + ',')
    return txt1.replace(',)', ')')


def dict2gauformat(dct: dict, gau_route=False):
    txt = ''
    if gau_route:
        txt = txt + '#p ' + dct['approach_basis']
        del dct['approach_basis']
    for upper_key, upper_val in dct.items():
        txt += (' ' + upper_key)
        if isinstance(upper_val, dict):
            txt += ('=(' + recursive_print(upper_val) + ')')
        elif bool(upper_val):
            txt += ('=' + upper_val)
    return txt
