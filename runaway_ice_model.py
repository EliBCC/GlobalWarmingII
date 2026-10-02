"""
Minimal version of script purely for autograder
"""
LRange = [1200, 1600]
albedoRange = [0.15, 0.65]  
iterations = 100
epsilon = 1
sigma = 5.67E-8          # W/m2 K4
max_iterations = 1000

m2 = -0.010000000000000016
b2 = 2.800000000000003

def tempFromAlbedo(L, albedo):
    return pow(L*(1-albedo)/(4*sigma*epsilon), 1/4)

def albedoFromTemp(temp):
    return m2*temp + b2

def albedoTempLoop(L, albedo, max_iterations = max_iterations):
    precision = 0.01
    old_albedo = albedo
    for _ in range(max_iterations):
        temp = tempFromAlbedo(L, albedo)
        albedo = albedoFromTemp(temp)
        # Hold albedo within the specified range
        albedo = max(min(albedo, albedoRange[1]), albedoRange[0])
        if (old_albedo - precision < albedo < old_albedo + precision):
            return temp, albedo
        old_albedo = albedo
    return temp, albedo

if __name__ == "__main__":
    L, albedo, nIters = input("").split()
    L, albedo, nIters = [ float(L), float(albedo), int(nIters) ]

    temp, albedo = albedoTempLoop(L, albedo, nIters)
    print(f"{temp} {albedo}")