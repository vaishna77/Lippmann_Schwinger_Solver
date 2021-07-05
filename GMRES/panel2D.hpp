class panel2D{
public:
	kernel_dtype integral;
	double xcenter;
	double ycenter;
	double xradius;
	double yradius;
	double* xNodes;
	double* yNodes;
	bool isLeaf;
	panel2D(double xcenter, double ycenter, double xradius, double yradius){
		this->xcenter	=	xcenter;
		this->ycenter	=	ycenter;
		this->xradius	=	xradius;
		this->yradius	=	yradius;
		this->isLeaf	=	false;
		this->integral	=	0.0+0.0*I;
	}
	~panel2D() {
		delete[] xNodes;
		delete[] yNodes;
	}
};
