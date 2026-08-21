import { isDirectRun, runCliScript } from './core/index.ts'
import { generateMetadataItems } from './generateProjects.ts'

/**
 * 课题为纯元信息模式：仅从 categories.yaml 读取课题元数据。
 * 与 projects 共用 generateProjects.ts 的通用实现 generateMetadataItems。
 */
function generateTopicsJson(locale = 'zh-CN') {
  generateMetadataItems('topics', locale)
}

if (isDirectRun(import.meta)) {
  runCliScript('topics.json', () => {
    generateTopicsJson('zh-CN')
    generateTopicsJson('en-US')
  })
}

export { generateTopicsJson }
