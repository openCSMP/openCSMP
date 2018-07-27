#ifndef TUTORIAL4_EXAMPLE_REVISITED_H
#define TUTORIAL4_EXAMPLE_REVISITED_H

#include "Example.h"
#include "Model.h"

namespace csmp {

	struct T4R_InputData {
		T4R_InputData() :set_(false) {};
		bool set_;
		std::string _model_name, _output_name;

	};

	class  Tutorial4_Example_Revisited : public Example
	{

	public:
		Tutorial4_Example_Revisited() {};
		Tutorial4_Example_Revisited(std::string, std::string);
		virtual void Run();
		virtual void Specifications();

	private:
		void scaleRegion(Model<2U>& sg, double64 scale_factor);
		void constructVelocityVector(Model<2U>& mdl);
		void assignFluxToPointSource(Model<2U>& mdl, const char* flux);
		T4R_InputData _inpdata;

	};

} // csmp

#endif // TUTORIAL4_EXAMPLE_REVISITED_H