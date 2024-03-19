#ifndef GUARD_SESSION_MANAGER
#define GUARD_SESSION_MANAGER

#include "Logger.hxx"
#include <memory>
#include <onnxruntime_cxx_api.h>
#include <string>
#include <unordered_map>

class OnnxSessionManager {
  public:
    OnnxSessionManager() {
        OrtLoggingLevel logging_level =
            ORT_LOGGING_LEVEL_WARNING; // ORT_LOGGING_LEVEL_VERBOSE

        env = Ort::Env(logging_level, "Default");
        session_options.SetInterOpNumThreads(1);
        session_options.SetIntraOpNumThreads(1);
    };
    Ort::Session *getSession(const std::string &modelPath) {
        // check if session already exists in the sessions map
        if (sessions_map.count(modelPath) == 0) {
            sessions_map[modelPath] = std::make_unique<Ort::Session>(
                env, modelPath.c_str(), session_options);
            Logger::get("OnnxSessionManager")
                ->info("Created session for model: {}", modelPath);
        } else {
            Logger::get("OnnxSessionManager")
                ->info("Session already exists for model: {}", modelPath);
        }
        return sessions_map[modelPath].get();
    };

  private:
    std::unordered_map<std::string, std::unique_ptr<Ort::Session>> sessions_map;
    Ort::Env env;
    Ort::SessionOptions session_options;
};

namespace onnxhelper {
void prepare_model(Ort::Session *session,
                   std::vector<int64_t> &input_node_dims,
                   std::vector<int64_t> &output_node_dims, int &num_input_nodes,
                   int &num_output_nodes);

std::vector<float> run_interference(Ort::Session* session,
                                    Ort::AllocatorWithDefaultOptions allocator,
                                    std::vector<float> &evt_input,
                                    std::vector<int64_t> input_node_dims,
                                    std::vector<int64_t> output_node_dims,
                                    const int num_input_nodes, const int num_output_nodes);

} // namespace onnxhelper

#endif /* GUARD_SESSION_MANAGER */