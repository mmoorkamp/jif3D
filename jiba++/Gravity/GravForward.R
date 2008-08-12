gravforward<-function(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos).C("CalcScalarForward",
  as.double(XSizes),
  as.integer(length(XSizes)),as.double(YSizes),as.integer(length(YSizes)),as.double(ZSizes),as.integer(length(ZSizes)),
  as.double(Densities),as.double(XMeasPos),as.double(YMeasPos),as.double(ZMeasPos),as.integer(length(XMeasPos)),
  GravAcceleration=double(length(XMeasPos)))
  
  
  
  gravtensorforward<-function(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos).C("CalcTensorForward",
  as.double(XSizes),
  as.integer(length(XSizes)),as.double(YSizes),as.integer(length(YSizes)),as.double(ZSizes),as.integer(length(ZSizes)),
  as.double(Densities),as.double(XMeasPos),as.double(YMeasPos),as.double(ZMeasPos),as.integer(length(XMeasPos)),
  GravAcceleration=double(length(XMeasPos)*9))

