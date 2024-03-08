import ROOT
import numpy as np

def SumW2Err (nll, mini):
    r2 = mini.save()
    cov2 = r2.covarianceMatrix()

    nll.applyWeightSquared(True)
    mini.hesse()
    r1 = mini.save()
    
    cov1 = r1.covarianceMatrix()
    mat1 = np.zeros([cov1.GetNrows(),cov1.GetNrows()])
    mat2 = np.zeros([cov2.GetNrows(),cov2.GetNrows()])
#     decomp = ROOT.Math.CholeskyDecompGenDim(cov1.GetNrows(), cov1)
#     decomp.Invert(cov1)
    for i in range(cov1.GetNrows()):
        for j in range(cov1.GetNrows()):
            mat1[i][j] = cov1[i][j]
    for i in range(cov2.GetNrows()):
        for j in range(cov2.GetNrows()):
            mat2[i][j] = cov2[i][j]
    incov = np.linalg.inv(mat1)
    print("V = ", mat2)
    print("c^-1 = ", incov)
    finalm = mat2.dot(incov.dot(mat2))
    print("Vc^-1V = ", finalm)
    err = np.zeros(cov2.GetNrows())
    for i in range(cov2.GetNrows()):
        err[i]=np.sqrt(finalm[i][i])

    # print ("NLL cov = ", r2.covQual())
    # print ("Weighted NLL cov = ", r1.covQual())
    nll.applyWeightSquared(False)
    return err

def Minimizer_NLL(nll, printLevel = -1, eps = 100, offSet = False, strategy = 0):
    mini = ROOT.RooMinimizer(nll)
    mini.setPrintLevel(printLevel)
    mini.setEps(eps)
    mini.setOffsetting(offSet)
    mini.setStrategy(strategy)
    mini.minimize("Minuit2","migrad")
    mini.hesse()
    r = mini.save()
    return r

