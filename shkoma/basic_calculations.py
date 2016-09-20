def mass_calculation(sequence):
    aa_array = list(sequence)
    db = open("../AA_table.csv", "r")
    mass = 0
    water_mass = 18.01528

    while True:
        line = db.readline()
        print(line)
        if line == '':
            break

        for line in db.readline():
            db_array = line.split(';')
            print(db_array)
            for letter in aa_array:
                if db_array[0] != "Amino_Acid":
                    if letter == db_array[2]:
                        mass += float(db_array[7])
                        print(mass + water_mass)
                        return mass


mass_calculation('ACDDEA')