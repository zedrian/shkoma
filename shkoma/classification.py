from pandas import DataFrame, Series, concat


def load_data(received_file_name, missed_file_name):
    received = DataFrame.from_csv(received_file_name)
    received['received'] = Series([1.0 for i in range(0, len(received))], index=received.index)
    missed = DataFrame.from_csv(missed_file_name)
    missed['received'] = Series([0.0 for i in range(0, len(missed))], index=missed.index)

    print('received peptides loaded: {0}'.format(len(received)))
    print('missed peptides loaded: {0}'.format(len(missed)))
    total = concat([received, missed])
    print('total peptides: {0}'.format(len(total)))

    return total
