/**
 * 「读 YAML → 兜底 → 归一化 → 写 JSON」生成器骨架，generateAbout / generateResources 共用。
 */
import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { readYaml, writeJsonFile } from './fs.ts'
import { logWriteSuccess } from './cli.ts'
import { logError, logWarn } from './log.ts'

export interface YamlJsonGeneratorOptions {
  locale?: string
  /** 源文件名（不带后缀，如 about / resources） */
  baseName: string
  /** YAML 缺失/解析失败时的兜底值 */
  fallback: unknown
  /** 归一化：raw YAML → 写盘数据 */
  normalize: (raw: unknown) => unknown
}

/** 生成 <baseName><suffix>.json：源按 locale 分层，输出恒为 content/ 下 */
export function generateJsonFromYaml(options: YamlJsonGeneratorOptions): void {
  const locale = options.locale ?? 'zh-CN'
  const suffix = localeSuffix(locale)
  const yamlPath = path.join(srcDirFor(locale), `${options.baseName}${suffix}.yaml`)
  const outputPath = path.join(contentDir, `${options.baseName}${suffix}.json`)

  let raw: unknown = options.fallback
  if (fs.existsSync(yamlPath)) {
    const yamlConfig = readYaml(yamlPath)
    if (yamlConfig !== null && yamlConfig !== undefined) {
      raw = yamlConfig
    } else {
      logWarn(`YAML 解析失败，使用回退值: ${yamlPath}`)
    }
  } else {
    logWarn(`源 YAML 不存在，使用回退值: ${yamlPath}`)
  }

  const result = options.normalize(raw)

  try {
    writeJsonFile(outputPath, result)
    logWriteSuccess(outputPath)
  } catch (error) {
    logError(`生成失败: ${outputPath}:`, error)
  }
}
