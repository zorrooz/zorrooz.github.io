/**
 * TOC 构建：从渲染容器读取标题层级生成树，并补齐缺失的标题 id。
 */
export interface TocNode {
  id: string
  text: string
  level: number
  children: TocNode[]
}

/** 读取标题文本（剔除锚点按钮与尾部 #） */
export function getHeadingText(h: Element): string {
  try {
    const clone = h.cloneNode(true) as Element
    clone.querySelectorAll('.heading-anchor')?.forEach((a) => a.remove())
    return (clone.textContent || '').replace(/\s*#\s*$/, '').trim()
  } catch {
    return (h.textContent || '').replace(/\s*#\s*$/, '').trim()
  }
}

/** 为无 id 的标题生成稳定 id（中文/字母/数字/连字符，冲突追加序号） */
export function ensureHeadingIds(headings: Element[]): void {
  headings.forEach((h) => {
    if (h.id) return
    let safeId = h.textContent
      .trim()
      .toLowerCase()
      .replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g, '')
      .replace(/\s+/g, '-')
    if (!safeId) safeId = `section-${Math.random().toString(36).substring(2, 9)}`

    let finalId = safeId
    let count = 1
    while (document.getElementById(finalId)) finalId = `${safeId}-${count++}`
    h.id = finalId
  })
}

/** 按 levels 构建两级目录树（顶级 + 子级） */
export function buildTocTree(headings: Element[], levels: number[]): TocNode[] {
  const levelSet = new Set(levels)
  const topLevel = Math.min(...levels)
  const secondLevel = topLevel + 1

  const tocList: TocNode[] = []
  let currentTop: TocNode | null = null
  for (const h of headings) {
    const level = parseInt(h.tagName.substring(1), 10)
    if (!levelSet.has(level)) continue
    const node: TocNode = { id: h.id, text: getHeadingText(h), level, children: [] }
    if (level === topLevel) {
      tocList.push(node)
      currentTop = node
    } else if (currentTop && level >= secondLevel) currentTop.children.push(node)
    else tocList.push(node)
  }
  return tocList
}
