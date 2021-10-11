intensity = []
absorbance = []
timepoints = []
temperature_error = []
# Add a list in every list, where every list is a channel
for i in range(0, 4):
    intensity.append([])
    absorbance.append([])
    timepoints.append([])
timepoints.append([])
intensity[0].append([0, 1, 2, 3])

print(intensity)
