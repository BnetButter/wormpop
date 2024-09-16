# utils.py
from sqlalchemy.sql import expression
from sqlalchemy.schema import DefaultClause

def get_column_default(column):
    if column.default is None:
        return None
    if isinstance(column.default, DefaultClause):
        if isinstance(column.default.arg, expression.Function):
            return None
        return column.default.arg
    return column.default.arg

def starve_from_l1_arrest(num_days):
    L = 98.09173354410315
    x0 = 17.484181571141384
    k = -0.47620124366404964

    def reverse_sigmoid(x, L, x0, k):
        return L / (1 + np.exp(-k * (x - x0)))

    return reverse_sigmoid(num_days, L, x0, k)

# Other helper functions can be added here
