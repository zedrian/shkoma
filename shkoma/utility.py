from sys import stdout


# simple function converting numpy.bytes_ to human-readable string
def b2str(bytes):
    return str(bytes)[2:-1]


# smart progress bar
def show_progress(label, width, percentage):
    progress = '['
    for i in range(0, width):
        if i / width < percentage:
            progress += '#'
        else:
            progress += ' '
    progress += '] {0:.1%}'.format(percentage)
    print('\r' + label + progress, end='')
    stdout.flush()
