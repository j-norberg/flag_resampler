#include <vector>
#include <cassert>
#include <cstring> // for memset
#include <string> // for debugging

#include <queue> // for job-queue
#include <mutex> // for job-queue
#include <chrono>

struct Job
{
	virtual ~Job(){};
	virtual void run() = 0;
	bool done = false; // set to true by worker
};

// just for main thread to give jobs to workers
// the workers never give them back, but mark the Job
class WorkQueue {
public:
	std::queue<Job*> _queue;
	std::mutex _mutex;

	void add(Job* j) {
//		puts("adding Job"); fflush(stdout);
		j->done = false;
		std::lock_guard<std::mutex> g(_mutex);
		_queue.push(j);
	}

	// return nullptr if empty
	Job* try_take() {
		Job* j = nullptr;
		std::lock_guard<std::mutex> g(_mutex);
		if (!_queue.empty())
		{
			j = _queue.front();
			_queue.pop();
		}
		return j;
	}

	// keep_running
	static void Worker(int worker_id, std::atomic_bool* keep_running, WorkQueue* q)
	{
//		bool debug = true;
		bool debug = false;

		// windows can not sleep for less than 14ms
		auto sleep_time = std::chrono::milliseconds(1);
		int consecutive_null = 0;
		while (*keep_running)
		{
			Job* j = q->try_take();
			if (j == nullptr)
			{
				if (consecutive_null < 10)
					std::this_thread::yield();
				else
					std::this_thread::sleep_for(sleep_time);

				if (debug)
				{
					printf("worker %d sleeping %d\n", worker_id, consecutive_null);
					fflush(stdout);
				}

				++consecutive_null;
				continue;
			}

			// restart count
			consecutive_null = 0;

//			printf("WORK! worker %d\n", worker_id);

			j->run();
			j->done = true;
		}
	}

};


