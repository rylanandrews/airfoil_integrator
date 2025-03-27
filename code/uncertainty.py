

def propagateRSS(func, values, uncertainties):
    uncertainty = 0

    for i in range(len(values)):
        # Approximate derivative of value i
        v1 = values
        v2 = values

        dx = values[i]*0.001
        v1[i] = v1[i] - dx
        v2[i] = v2[i] + dx

        df = (func(v2) - func(v1)) / (2*dx)

        uncertainty += (uncertainties[i] * df)^2

    return uncertainty

