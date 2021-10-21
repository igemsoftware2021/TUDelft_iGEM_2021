from file_helpers import read_absorbance_csv, write_absorbance_csv

time = [[0.1, 0.2, 0.35, 0.56], [0.1, 0.35, 0.56], [0.1, 0.25, 0.35, 0.56, 0.2]]
absorbance = [[0.2, 0.3, 0.43, 8], [0.2, 0.43, 8], [0.2, 0.35, 0.43, 8, 6.8]]

write_absorbance_csv(time, absorbance)
#time_r, absorbance_r = read_absorbance_csv(path="dbsorbance-1631875086.csv")
# print(time_r)
# print(absorbance_r)
