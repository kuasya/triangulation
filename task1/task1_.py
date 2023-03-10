import sys
import collections
stdout = sys.stdout 
data = {}

# сумма по всем каналам
for line in open('krf20090408_71203_2_S1.txt', 'r'):
    if line == 'PLAINFO\n': break
    row = line.split()
    interval = (round(float(row[0]), 3), round(float(row[1]), 3))
    count = sum(map(int, row[2:]))
    data[interval] = count

out_file = open('out_single.txt', 'w', encoding = 'utf-8')
sys.stdout = out_file
for interval, count in data.items():
    print('%.3f   %.3f   %s' % (interval[0], interval[1], count))
sys.stdout = stdout

# объединение интервалов до 0.064
intervals = list(data.keys())
union_intervals = {} 
deltas = []
# delete - 28 строк
passing = 28
passed = 0
interval_start = 0
interval_length = 0
interval_max  = 0.064
interval_value = 0
for interval in intervals:
    if passed < passing:
        passed += 1
        continue
    if not interval_start: interval_start = interval[0]
    delta = round(interval[1] - interval[0], 3)
    interval_length = round(interval_length + delta, 3)
    interval_value  += data[interval]
    if interval_length + delta > interval_max + 0.001:
        interval_end = interval_start + interval_length
        union_intervals[(interval_start, interval_end)] = interval_value
        deltas.append(interval_length)
        interval_start  = interval_end
        interval_length = 0
        interval_value  = 0

counter = collections.Counter(deltas)
print(counter)

out_file = open('krf_64.txt', 'w', encoding = 'utf-8')
sys.stdout = out_file
for interval, count in union_intervals.items():
    print('%.3f   %.3f   %s' % (interval[0], interval[1], count))

sys.stdout = stdout

