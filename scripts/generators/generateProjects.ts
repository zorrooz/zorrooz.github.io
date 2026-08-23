import { isDirectRun, runCliScript } from '../lib/cli.ts'
import { generateMetadataItems } from '../lib/metadata.ts'

/**
 * 项目为纯元信息模式：不从 md 扫描文章，仅从 categories.yaml 读取项目元数据。
 * 实现已下沉 scripts/lib/metadata.ts（与 topics 共用）。
 */
function generateProjectsJson(locale = 'zh-CN') {
  generateMetadataItems('projects', locale)
}

if (isDirectRun(import.meta)) {
  runCliScript('projects.json', () => {
    generateProjectsJson('zh-CN')
    generateProjectsJson('en-US')
  })
}

export { generateProjectsJson }
