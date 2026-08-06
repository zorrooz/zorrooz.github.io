/**
 * 领域类型：与数据分支 content/*.json 产物一一对应（见 scripts/generators/）。
 * 运行时消费统一从这里取类型，不再在各 view 内散落 any。
 */

export interface Post {
  id: number
  no: number
  title: string
  date: string
  category: [string, string]
  tags: string[]
  preview: string
  wordCount: number
}

export interface Note {
  title: string
  date: string
  tags: string[]
  draft: boolean
  description: string
  relativePath: string
  wordCount: number
  tagCount: number
}

export interface Tag {
  name: string
  count: number
}

export interface CategoryArticle {
  title: string
  articleUrl: string
  wordCount: number
  tags: string[]
}

export interface CategoryStats {
  postsCount: number
  totalWords: number
  latestDate: string
}

export interface CategorySubcategory {
  key: string
  title: string
  articles: CategoryArticle[]
  stats: CategoryStats
}

export interface CategoryItem {
  name: string
  title: string
  desc: string
  tags: string[]
  stats: CategoryStats
  root?: string
  categories: CategorySubcategory[]
  /** notes 分类直属文章（当前产物未输出，保留防御） */
  articles?: CategoryArticle[]
  /** projects / topics 元数据 */
  github?: string
  url?: string
  doi?: string
  status?: string
  language?: string
  stars?: number
  license?: string
  version?: string
  journal?: string
  year?: number
  authors?: string[]
}

export interface CategorySection {
  title: string
  items: CategoryItem[]
}

export type CategoryData = CategorySection[]

export interface ResourceItem {
  name: string
  url: string
  desc: string
}

export interface ResourceNode {
  title: string
  items?: ResourceItem[]
  children?: ResourceNode[]
}

export interface AboutExperience {
  year: string
  title: string
  desc: string
}

export interface AboutSectionItem {
  name: string
  desc: string
}

export interface AboutSection {
  title: string
  items: AboutSectionItem[]
}

export interface AboutContact {
  label: string
  value: string
  link: string
  icon: string
}

export interface AboutData {
  introduction: string
  experience: AboutExperience[]
  section: AboutSection[]
  contacts: AboutContact[]
}
