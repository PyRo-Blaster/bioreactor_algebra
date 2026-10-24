import numpy as np, itertools
C = np.array([[0,1,1,1,1,1],[1,0,1,-1,-1,1],[1,1,0,1,-1,-1],[1,-1,1,0,1,-1],[1,-1,-1,1,0,1],[1,1,-1,-1,1,0]],float)
print("C^T C = 5I:", np.allclose(C.T@C, 5*np.eye(6)))
D = np.vstack([C,-C,np.zeros((1,6))])          # 13-run DSD, 6 factors
def phi(x):
    q=[x[i]*x[j] for i,j in itertools.combinations(range(6),2)]
    return np.concatenate([[1],x,x**2,q])
X = np.array([phi(r) for r in D]); n,P = X.shape
print("runs, terms:", n, P, " rank(X):", np.linalg.matrix_rank(X), " null-space dim:", P-np.linalg.matrix_rank(X))
# correlations
Xc = X[:,1:]-X[:,1:].mean(0); sd = Xc.std(0); R = (Xc/sd).T@(Xc/sd)/n
lin=range(0,6); quad=range(6,12); inter=range(12,27)
print("max |r| linear vs 2nd-order:", np.abs(R[np.ix_(lin,list(quad)+list(inter))]).max().round(3))
print("r among quadratic columns:", np.unique(R[np.ix_(quad,quad)][np.triu_indices(6,1)].round(3)))
print("max |r| quad vs 2FI:", np.abs(R[np.ix_(quad,inter)]).max().round(3), "  max |r| among 2FI:", np.abs(R[np.ix_(inter,inter)][np.triu_indices(15,1)]).max().round(3))
# estimability of corner predictions (is phi(x) in row space of X?)
Q,_ = np.linalg.qr(X.T)            # basis of row space
Q = Q[:, :np.linalg.matrix_rank(X)]
def unexplained(x):
    f=phi(np.array(x,float)); r=f-Q@(Q.T@f); return np.linalg.norm(r)/np.linalg.norm(f)
print("share of phi(x) outside estimable space:")
print("  centre          ", round(unexplained([0]*6),3))
print("  a DSD run       ", round(unexplained(D[0]),3))
print("  corner (+1..+1) ", round(unexplained([1]*6),3))
print("  corner (-1..-1) ", round(unexplained([-1]*6),3))
print("  mixed corner    ", round(unexplained([1,-1,1,-1,1,-1]),3))
# pass-probability gate if criterion = inside own 95% predictive interval
for m in (4,6,8): print(f"m={m} responses x 2 vessels, independent: P(all inside 95% PI) = {0.95**(2*m):.2f}")
# corrected per-interval level for a joint 95%
for m in (6,): print("Bonferroni per-interval level for joint 0.95, m=6, 2 vessels:", round(1-0.05/(2*m),4))
# superficial gas velocity ratio at constant vvm, geometric similarity 7 -> 15 L
print("v_s ratio 15L/7L at equal vvm (geom. similar):", round((15/7)**(1/3),3))
