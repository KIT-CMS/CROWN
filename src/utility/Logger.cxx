#include <utility/Logger.hxx>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <vector>

std::shared_ptr<spdlog::logger> Logger::get(std::string name) {
    if (getInstance()._loggers.count(name) == 0) {
        std::vector<spdlog::sink_ptr> sinkVector;
        sinkVector.push_back(
            std::make_shared<spdlog::sinks::stdout_color_sink_st>());

        if (getInstance()._fileName)
            sinkVector.push_back(
                std::make_shared<spdlog::sinks::basic_file_sink_st>(
                    *getInstance()._fileName));

        auto newLogger = std::make_shared<spdlog::logger>(
            name, begin(sinkVector), end(sinkVector));
        newLogger->set_level(convertLevelToSpdlog(getInstance()._level));
        getInstance()._loggers[name] = newLogger;
    }

    return getInstance()._loggers[name];
}

void Logger::setLevel(LogLevel level) {
    getInstance()._level = level;
    spdlog::set_level(convertLevelToSpdlog(level));

    for (auto &[key, logger] : getInstance()._loggers)
        logger->set_level(convertLevelToSpdlog(level));
}

void Logger::enableFileLogging(std::string filename) {
    getInstance()._fileName = std::make_unique<std::string>(filename);
    for (auto &[key, logger] : getInstance()._loggers) {
        if (logger->sinks().size() < 2)
            logger->sinks().push_back(
                std::make_shared<spdlog::sinks::basic_file_sink_st>(
                    *getInstance()._fileName));
    }
}

Logger::~Logger() {
    _loggers.clear();
    _fileName.reset();
}

Logger &Logger::getInstance() {
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
