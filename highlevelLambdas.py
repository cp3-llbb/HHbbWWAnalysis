from bamboo import treefunctions as op

# 4-Momentum association #
ll_p4 = lambda l1,l2 : l1.p4+l2.p4
lljj_p4 = lambda l1,l2,j1,j2 : l1.p4+l2.p4+j1.p4+j2.p4

# Dilep-Met variables #
DilepMET_deltaPhi = lambda l1,l2,met : ll_p4(l1,l2).Phi()-met.phi
DilepMET_Pt = lambda l1,l2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+ll_p4(l1,l2).Px(),2)+op.pow(met.pt*op.sin(met.phi)+ll_p4(l1,l2).Py(),2))

# Transverse mass #
MT_ll = lambda l1,l2,met : op.sqrt(2*ll_p4(l1,l2).Pt()*met.pt*(1-op.cos(ll_p4(l1,l2).Phi()-met.phi)))
MT_lljj = lambda l1,l2,j1,j2,met : op.sqrt(2*lljj_p4(l1,l2,j1,j2).Pt()*met.pt*(1-op.cos(lljj_p4(l1,l2,j1,j2).Phi()-met.phi)))

# dilep + dijet #
#M_lljj = lambda l1,l2,j1,j2 : op.invariant_mass(ll_p4(l1,l2),jj_p4(j1,j2))
M_lljj = lambda l1,l2,j1,j2 : op.invariant_mass(lljj_p4(l1,l2,j1,j2))
MinDR_lj = lambda l1,l2,j1,j2 : op.min(op.min(op.deltaR(l1.p4,j1.p4),op.deltaR(l1.p4,j2.p4)),
                                       op.min(op.deltaR(l2.p4,j1.p4),op.deltaR(l2.p4,j2.p4)))

# Higgs related variables #
HT2 = lambda l1,l2,j1,j2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l1.p4.Px()+l2.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l1.p4.Py()+l2.p4.Py(),2)) + op.abs((j1.p4+j2.p4).Pt())
HT2R = lambda l1,l2,j1,j2,met : HT2(met,l1,l2,j1,j2)/(met.pt+l1.p4.Pt()+l2.p4.Pt()+j1.p4.Pt()+j2.p4.Pt())
