/**
 * 文本清洗共享实现（自 scripts/generators/core 与 generateSearchIndex 收敛）。
 * markdownToPlain 与 stripHtmlToText 遵循同一风格：剥除标记/标签/实体 → 折叠空白 → trim。
 */

/** 粗略去除 markdown 标记 → 纯文本，用于字数统计 */
export function markdownToPlain(text: string): string {
  // 先剥代码块与行内代码
  let t = text.replace(/```[\s\S]*?```/g, ' ')
  t = t.replace(/`[^`]*`/g, ' ')
  // 图像/链接/HTML 标签
  t = t.replace(/!\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/<[^>]+>/g, ' ')
  // 标题、加粗/斜体/引用/列表/表格标记与参考链接
  t = t.replace(/^#{1,6}\s+/gm, ' ')
  t = t.replace(/[*_~`>#|-]{1,}/g, ' ')
  t = t.replace(/^\s*\[[^\]]+]:\s+\S+.*$/gm, ' ')
  // 折叠空白
  t = t.replace(/\s+/g, ' ').trim()
  return t
}

export function countWordsSmart(text: string): number {
  const cjk = (text.match(/[\u4E00-\u9FFF\u3400-\u4DBF]/g) || []).length
  const words = (text.match(/[A-Za-z0-9]+(?:'[A-Za-z0-9]+)?/g) || []).length
  return cjk + words
}

/** HTML → 纯文本：剥除 style/script 与标签、常见实体解码后折叠空白（搜索索引用） */
export function stripHtmlToText(html: string): string {
  return html
    .replace(/<style[\s\S]*?<\/style>/gi, ' ')
    .replace(/<script[\s\S]*?<\/script>/gi, ' ')
    .replace(/<[^>]+>/g, ' ')
    .replace(/&amp;/gi, '&')
    .replace(/&lt;/gi, '<')
    .replace(/&gt;/gi, '>')
    .replace(/&quot;/gi, '"')
    .replace(/&#39;|&apos;/gi, "'")
    .replace(/&nbsp;/gi, ' ')
    .replace(/\s+/g, ' ')
    .trim()
}
