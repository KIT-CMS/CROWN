#ifndef GUARDLOGGER_H
#define GUARDLOGGER_H

#include <map>
#include <spdlog/fmt/ostr.h> // for formatting of RVecs
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

class Logger {
  public:
    static std::shared_ptr<spdlog::logger> get(std::string name) {
        if (getInstance()._loggers.count(name) == 0) {
            std::vector<spdlog::sink_ptr> sinkVector;
            sinkVector.push_back(
                std::make_shared<spdlog::sinks::stdout_color_sink_st>());

            // check if file logging is enabled
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
    enum class LogLevel { DEBUG, INFO, WARN, ERR, CRITICAL, OFF };
    static void setLevel(LogLevel level) {
        getInstance()._level = level;

        // set level globally (probably superfluous..)
        spdlog::set_level(convertLevelToSpdlog(level));

        // set level for all active loggers
        for (auto &[key, logger] : getInstance()._loggers)
            logger->set_level(convertLevelToSpdlog(level));
    }
    static void enableFileLogging(std::string filename) {
        getInstance()._fileName = std::make_unique<std::string>(filename);
        for (auto &[key, logger] : getInstance()._loggers) {
            // if there is less than two sinks, add a file sink
            if (logger->sinks().size() < 2)
                logger->sinks().push_back(
                    std::make_shared<spdlog::sinks::basic_file_sink_st>(
                        *getInstance()._fileName));
        }
    }

  private:
    static Logger &getInstance() {
        static Logger instance;
        return instance;
    }
    LogLevel _level{LogLevel::INFO};
    static spdlog::level::level_enum convertLevelToSpdlog(LogLevel level) {
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
    std::unique_ptr<std::string> _fileName{};
    std::map<std::string, std::shared_ptr<spdlog::logger>> _loggers;
};

#endif /* GUARDLOGGER_H */
