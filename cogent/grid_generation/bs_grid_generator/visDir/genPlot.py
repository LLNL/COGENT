DeleteAllPlots()
DefineVectorExpression("disp", "{component_0,component_1}-coords(Mesh)")
path = "Results/ghostCellExtrapBS/diagnostic/"
filenames = [ "16.aligned.final.lcore.2d.hdf5", "16.aligned.final.lcsol.2d.hdf5", "16.aligned.final.lpf.2d.hdf5", "16.aligned.final.lsol.2d.hdf5", "16.aligned.final.mcore.2d.hdf5", "16.aligned.final.mcsol.2d.hdf5", "16.aligned.final.rcore.2d.hdf5", "16.aligned.final.rcsol.2d.hdf5", "16.aligned.final.rpf.2d.hdf5", "16.aligned.final.rsol.2d.hdf5" ]
for filename in filenames:
   OpenDatabase(path + filename)
   AddPlot("Mesh", "Mesh")
   x = GetPlotOptions()
   x.legendFlag = 0
   SetPlotOptions(x)
DrawPlots()
AddOperator("AMRStitchCell", 1)
AddOperator("DeferExpression", 1)
DeferExpressionAtts = DeferExpressionAttributes()
DeferExpressionAtts.exprs = ("disp")
SetOperatorOptions(DeferExpressionAtts, 1, 1)
AddOperator("Displace", 1)
DisplaceAtts = DisplaceAttributes()
DisplaceAtts.factor = 1
DisplaceAtts.variable = "disp"
SetOperatorOptions(DisplaceAtts, 1, 1)
DrawPlots()
ResetView()

