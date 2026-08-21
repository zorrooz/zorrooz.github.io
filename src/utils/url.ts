/** DOI 正则（全站唯一一份） */
const DOI_PATTERN = /^10\.\d{4,9}\//

/** 规范化外部链接：空值返回 ''，http(s) 原样返回，其余补 https:// 前缀并去开头斜杠 */
export function normalizeUrl(value: string | undefined): string {
  if (!value || !value.trim()) return ''
  if (/^https?:\/\//i.test(value)) return value
  return 'https://' + value.replace(/^\/+/, '')
}

/** 规范化 DOI：空值返回 ''，http(s) 原样返回，裸 DOI 转 doi.org 链接，其余补 https:// */
export function normalizeDoi(value: string | undefined): string {
  if (!value || !value.trim()) return ''
  if (/^https?:\/\//i.test(value)) return value
  if (DOI_PATTERN.test(value)) return 'https://doi.org/' + value
  return 'https://' + value.replace(/^\/+/, '')
}

/** 展示用 URL：去掉协议前缀与末尾斜杠 */
export function displayUrl(url: string): string {
  if (!url) return ''
  return url.replace(/^https?:\/\//, '').replace(/\/$/, '')
}

/** 判断 URL 是否为 DOI 链接（含 doi.org 或裸 DOI） */
export function isDoi(url: string): boolean {
  return !!url && (url.includes('doi.org') || DOI_PATTERN.test(url))
}
