#include <utility/Logger.hxx>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <vector>
#include <string_view>


static const std::shared_ptr<spdlog::logger>& Logger::get(std::string_view name) {
    // per-thread cache keeps the hot path lock-free
    thread_local std::unordered_map<std::string_view,
                                    const std::shared_ptr<spdlog::logger> *>
        cache;
    auto cached = cache.find(name);
    if (cached != cache.end())
        return *cached->second;

    Logger &instance = getInstance();
    const std::string key(name);
    std::unique_lock<std::shared_mutex> lock(instance._mutex);
    auto it = instance._loggers.find(key);
    if (it == instance._loggers.end()) {
        std::vector<spdlog::sink_ptr> sinkVector;
        sinkVector.push_back(
            std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        if (instance._fileName)
            sinkVector.push_back(
                std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                    *instance._fileName));
        auto newLogger = std::make_shared<spdlog::logger>(
            key, begin(sinkVector), end(sinkVector));
        newLogger->set_level(convertLevelToSpdlog(instance._level));
        it = instance._loggers.emplace(key, std::move(newLogger)).first;
    }
    cache.emplace(std::string_view(it->first), &it->second);
    return it->second;
}


static void Logger::setLevel(LogLevel level) {
    std::unique_lock<std::shared_mutex> writeLock(getInstance()._mutex);
    getInstance()._level = level;

    // set level globally (probably superfluous..)
    spdlog::set_level(convertLevelToSpdlog(level));

    // set level for all active loggers
    for (auto &[key, logger] : getInstance()._loggers)
        logger->set_level(convertLevelToSpdlog(level));
}


static void enableFileLogging(std::string filename) {
    std::unique_lock<std::shared_mutex> writeLock(getInstance()._mutex);
    getInstance()._fileName = std::make_unique<std::string>(filename);
    for (auto &[key, logger] : getInstance()._loggers) {
        // if there is less than two sinks, add a file sink
        if (logger->sinks().size() < 2)
            logger->sinks().push_back(
                std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                    *getInstance()._fileName));
    }
}


Logger::~Logger() {
    _loggers.clear();
    _fileName.reset();
}


static Logger &Logger::getInstance() {
    static Logger instance;
    return instance;
}


spdlog::level::level_enum Logger::convertLevelToSpdlog(LogLevel level) {
    switch (level) {
    case LogLevel::DEBUG:
        return spdlog::level::debug;
    case LogLevel::INFO:
        return spdlog::level::info;
    case LogLevel::WARN:
        return spdlog::level::warn;
    case LogLevel::ERR:
        return spdlog::level::err;
    case LogLevel::CRITICAL:
        return spdlog::level::critical;
    case LogLevel::OFF:
        return spdlog::level::off;
    default:
        return spdlog::level::info;
    }
}
