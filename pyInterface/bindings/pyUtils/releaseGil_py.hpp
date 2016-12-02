#ifndef RELEASEGIL_PY_H
#define RELEASEGIL_PY_H


namespace rpwa {


	namespace py {


		class releaseGil {

		public:

			releaseGil()
			{
				state = PyEval_SaveThread();
			}

			~releaseGil()
			{
				PyEval_RestoreThread(state);
			}

		private:

			PyThreadState* state;

		};


	}


}


#endif // RELEASEGILE_PY_H
