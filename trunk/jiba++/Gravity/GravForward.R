AllocateModel<-function(cachescalar,cachetensor).C("AllocateModel",as.integer(cachescalar),as.integer(cachetensor))

GravScalarForward<-function(XSizes,YSizes,ZSizes,Densities,BG_Densities,BG_Thicknesses,XMeasPos,YMeasPos,ZMeasPos).C("CalcScalarForward",
  as.double(XSizes),
  as.integer(length(XSizes)),as.double(YSizes),as.integer(length(YSizes)),as.double(ZSizes),as.integer(length(ZSizes)),
  as.double(Densities),as.double(BG_Densities),as.double(BG_Thicknesses),as.integer(length(BG_Thicknesses)),
  as.double(XMeasPos),as.double(YMeasPos),as.double(ZMeasPos),as.integer(length(XMeasPos)),
  GravAcceleration=double(length(XMeasPos)))
  
  
  
GravTensorForward<-function(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos).C("CalcTensorForward",
  as.double(XSizes),
  as.integer(length(XSizes)),as.double(YSizes),as.integer(length(YSizes)),as.double(ZSizes),as.integer(length(ZSizes)),
  as.double(Densities),as.double(BG_Densities),as.double(BG_Thicknesses),as.integer(length(BG_Thicknesses)),
  as.double(XMeasPos),as.double(YMeasPos),as.double(ZMeasPos),as.integer(length(XMeasPos)),
  GravAcceleration=double(length(XMeasPos)*9))

