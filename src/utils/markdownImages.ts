/**
 * 文章图片资源解析：@data 别名（zh 源 content-src / 英文镜像 cache/en）静态资源 glob，
 * 将渲染 HTML 里相对路径的 <img> 重写为构建产物 URL，兼容 Obsidian「最短路径」引用。
 */

const assetModules = import.meta.glob(
  '@data/{content-src,cache/en}/**/*.{png,jpg,jpeg,gif,svg,webp}',
  { query: '?url', import: 'default', eager: true },
)

/** 相对文章目录解析资源；找不到时返回原路径 */
const resolveAssetUrl = (relPath: string, articleDir: string): string => {
  if (/^(https?:)?\/\//i.test(relPath) || relPath.startsWith('/')) return relPath

  const findAsset = (normalized: string, keys: string[]): string | null => {
    for (const key of keys) {
      const matched = Object.keys(assetModules).find((k) => k.endsWith(`/${normalized}`) || k === key)
      if (matched) return assetModules[matched] as string
    }
    return null
  }

  const parts = (articleDir + '/' + relPath).split('/').filter((p) => p && p !== '.')
  const stack: string[] = []
  parts.forEach((p) => (p === '..' ? stack.pop() : stack.push(p)))
  const normalized = stack.join('/')
  const byDir = findAsset(normalized, [`@data/content-src/${normalized}`, normalized])
  if (byDir) return byDir

  /** Obsidian 最短路径：以 vault 根为基准（如 assets/fig-1.png） */
  const rootRef = relPath.replace(/^\.\//, '')
  if (rootRef.startsWith('assets/')) {
    const byRoot = findAsset(rootRef, [`@data/content-src/${rootRef}`, rootRef])
    if (byRoot) return byRoot
  }
  return relPath
}

/** 重写 HTML 中全部 <img src>；解析异常时返回原文并告警 */
export function rewriteImageLinks(html: string, articlePath: string): string {
  try {
    const articleDir = articlePath
      .replace(/^[./]*/, '')
      .replace(/\.md$/, '')
      .split('/')
      .slice(0, -1)
      .join('/')

    return html.replace(
      /<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,
      (_m, pre, src, post) => `<img ${pre}src="${resolveAssetUrl(src.trim(), articleDir)}"${post}>`,
    )
  } catch (e) {
    console.warn('rewriteImageLinks failed', e)
    return html
  }
}
