"""Aggregate FIA species-group subtotals into 11 fall-foliage categories."""
import argparse
import csv
import math
from pathlib import Path

GROUPS = {
    'Evergreen': ['Longleaf and slash pines', 'Loblolly and shortleaf pines', 'Other yellow pines', 'Eastern white and red pines', 'Jack pine', 'Spruce and balsam fir', 'Eastern hemlock', 'Other eastern softwoods', 'Urban-specific softwoods'],
    'Soft maple': ['Soft maple'],
    'Hard maple': ['Hard maple'],
    'Oak': ['Select white oaks', 'Select red oaks', 'Other white oaks', 'Other red oaks'],
    'Beech': ['Beech'],
    'Sweetgum': ['Sweetgum'],
    'Yellow-poplar': ['Yellow-poplar'],
    'Tupelo': ['Tupelo and blackgum'],
    'Golden hardwood': ['Hickory', 'Yellow birch', 'Cottonwood and aspen (East)', 'Basswood', 'Black walnut'],
    'Ash': ['Ash'],
    'Miscellaneous deciduous': ['Cypress', 'Other eastern soft hardwoods', 'Other eastern hard hardwoods', 'Eastern noncommercial hardwoods', 'Urban-specific hardwoods'],
}
MODEL_GROUPS = {name.casefold(): category for category, names in GROUPS.items() for name in names}

def distribution(path):
    counts = dict.fromkeys(GROUPS, 0.0)
    seen = set()
    total = None
    with open(path, encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        required = {'SPECIES_GROUP', 'ALL_LIVE_STOCKING', 'ESTIMATE'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('CSV must contain SPECIES_GROUP, ALL_LIVE_STOCKING, ESTIMATE')
        for row in reader:
            if row['ALL_LIVE_STOCKING'].strip().casefold() != 'total':
                continue
            name = row['SPECIES_GROUP'].strip()
            key = name.casefold()
            if key in seen:
                raise ValueError(f'Duplicate subtotal: {name}')
            seen.add(key)
            value = float(row['ESTIMATE'])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'Invalid estimate for {name}: {value}')
            if key == 'total':
                total = value
            elif key not in MODEL_GROUPS:
                raise ValueError(f'Unclassified FIA group: {name}. Add it to GROUPS before proceeding.')
            else:
                counts[MODEL_GROUPS[key]] += value
    if total is None or total <= 0:
        raise ValueError('Missing or nonpositive grand total')
    if not math.isclose(sum(counts.values()), total, rel_tol=1e-8, abs_tol=0.01):
        raise ValueError(f'Group subtotals ({sum(counts.values())}) do not match grand total ({total})')
    return [(name, count, 100 * count / total) for name, count in counts.items()]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv_file', type=Path, help='Input city FIA CSV')
    parser.add_argument('-o', '--output', type=Path, help='Output CSV (default: beside input)')
    args = parser.parse_args()
    output = args.output or args.csv_file.with_name(args.csv_file.stem + '_model_distribution.csv')
    try:
        if output.resolve() == args.csv_file.resolve():
            raise ValueError('Output must not overwrite the input CSV')
        rows = distribution(args.csv_file)
        with output.open('w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['MODEL_CATEGORY', 'ESTIMATED_TREES', 'PERCENTAGE'])
            writer.writerows(rows)
        print(f'{"Category":30} {"Estimated trees":>18} {"Percent":>10}')
        for name, count, percent in rows:
            print(f'{name:30} {count:18,.0f} {percent:9.2f}%')
        print(f'{"Total":30} {sum(r[1] for r in rows):18,.0f} {sum(r[2] for r in rows):9.2f}%')
        print(f'\nSaved: {output}')
    except (OSError, ValueError) as e:
        parser.exit(1, f'Error: {e}\n')

if __name__ == '__main__':
    main()
