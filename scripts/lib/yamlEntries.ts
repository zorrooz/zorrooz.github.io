import { safeArray } from '../generators/core/index.ts'

/**
 * categories.yaml 的 projects/topics 条目统一归一化。
 * 收敛 generateProjects / generateTopics / generateCategories 三处重复的
 * name/title/desc/date/tags 等 typeof 防御写法；字段集 = 三处现有字段的并集。
 *
 * 语义约定（保证与原三处实现逐字节等价）：
 * - name 非字符串 → ''，无效条目（空 name）由调用方丢弃
 *   （generateProjects/generateTopics 原以 `typeof name !== 'string' || !name` 置 null，
 *   generateCategories 原置 '' 后以 `!config.name` 丢弃，二者收敛于同一结果）
 * - title 非字符串 → 回退为 name（generateProjects/generateTopics 原逻辑；
 *   generateCategories 以 `title || name` 后处理，结果一致）
 * - tags/authors 用数组过滤（safeArray + 字符串元素过滤），与三处原实现一致
 */
export type NormalizedProjectTopicEntry = {
  name: string
  title: string
  desc: string
  date: string
  tags: string[]
  github: string
  doi: string
  url: string
  status: string
  language: string
  stars: number
  license: string
  version: string
  journal: string
  year: number
  authors: string[]
}

export function normalizeProjectTopicEntry(raw: unknown): NormalizedProjectTopicEntry {
  const def = (raw ?? {}) as Record<string, unknown>
  return {
    name: typeof def?.name === 'string' ? def.name : '',
    title:
      typeof def?.title === 'string' ? def.title : typeof def?.name === 'string' ? def.name : '',
    desc: typeof def?.desc === 'string' ? def.desc : '',
    date: typeof def?.date === 'string' ? def.date : '',
    tags: safeArray(def?.tags).filter((t) => typeof t === 'string'),
    github: typeof def?.github === 'string' ? def.github : '',
    doi: typeof def?.doi === 'string' ? def.doi : '',
    url: typeof def?.url === 'string' ? def.url : '',
    status: typeof def?.status === 'string' ? def.status : '',
    language: typeof def?.language === 'string' ? def.language : '',
    stars: Number(def?.stars) || 0,
    license: typeof def?.license === 'string' ? def.license : '',
    version: typeof def?.version === 'string' ? def.version : '',
    journal: typeof def?.journal === 'string' ? def.journal : '',
    year: Number(def?.year) || 0,
    authors: safeArray(def?.authors).filter((t) => typeof t === 'string'),
  }
}
