<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="tree-group">
      <div class="tree-label">{{ category.name }}</div>
      <ul class="article-list">
        <TreeNode
          v-for="dir in category.children"
          :key="dir.id"
          :dir="dir"
          :depth="1"
          :collapsed-ids="collapsedIds"
          :is-active="isActive"
          :to-article="toArticle"
          @toggle="toggleDir"
        />

        <li v-for="file in category.files" :key="file.path" class="article-item">
          <router-link
            :to="toArticle(file.path)"
            class="tree-item tree-item--l1"
            :class="{ 'tree-item--active': isActive(file.path) }"
          >
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import { useRoute } from 'vue-router'
import { ARTICLE_ROUTE_PREFIX } from '@/config'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { loadCategories } from '@/utils/contentLoader'
import {
  articlePathFromUrl,
  joinRoutePathParam,
  normalizeArticleKey,
  toArticle,
} from '@/utils/navigation'
import TreeNode, { type TreeDir, type TreeFile } from './TreeNode.vue'
import type { CategoryItem } from '@/types'

interface TreeCategory {
  name: string
  files: TreeFile[]
  children: TreeDir[]
}

const route = useRoute()

const navigationTree = ref<TreeCategory[]>([])
const currentPath = ref('')
const collapsedIds = ref<string[]>([])
const { data: categoryData } = useLocalizedContent(() => loadCategories(), [])

function isActive(path: string) {
  return normalizeArticleKey(currentPath.value) === normalizeArticleKey(path)
}

/** 默认全部展开；仅记录用户手动折叠的目录 */
function isExpanded(dir: TreeDir) {
  return !collapsedIds.value.includes(dir.id)
}

function toggleDir(dir: TreeDir) {
  collapsedIds.value = isExpanded(dir)
    ? [...collapsedIds.value, dir.id]
    : collapsedIds.value.filter((id) => id !== dir.id)
}

function buildTree() {
  const path = currentPath.value
  if (!path) {
    navigationTree.value = []
    return
  }
  const segs = path.split('/').filter(Boolean)
  if (segs.length < 2) {
    navigationTree.value = []
    return
  }
  const type = segs[0]
  const group = segs[1]

  let targetItem: CategoryItem | null = null
  outer: for (const section of categoryData.value) {
    for (const item of section.items) {
      if (item.name !== group) continue

      const inArticle = (a: { articleUrl: string }) =>
        a.articleUrl.includes(`${ARTICLE_ROUTE_PREFIX}/${type}/`)
      const hasTypeMatch =
        item.articles?.some(inArticle) ||
        item.categories.some((cat) =>
          cat.articles.some((a) =>
            a.articleUrl.includes(`${ARTICLE_ROUTE_PREFIX}/${type}/${group}/`),
          ),
        )

      if (hasTypeMatch) {
        targetItem = item
        break outer
      }
    }
  }
  if (!targetItem) {
    navigationTree.value = []
    return
  }

  const rootFiles: TreeFile[] = []
  const children: TreeDir[] = []

  const toFile = (title: string, articleUrl: string): TreeFile => ({
    title,
    path: articlePathFromUrl(articleUrl),
  })

  targetItem.articles?.forEach((a) => {
    if (a.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl))
  })

  targetItem.categories.forEach((cat) => {
    const files: TreeFile[] = []
    cat.articles.forEach((a) => {
      if (a.articleUrl) files.push(toFile(a.title, a.articleUrl))
    })
    if (files.length) {
      children.push({
        id: `${targetItem!.name}/${cat.key || cat.title}`,
        name: cat.title || cat.key,
        files,
      })
    }
  })

  const tree: TreeCategory[] = [
    {
      name: targetItem.title || targetItem.name || group,
      files: rootFiles,
      children,
    },
  ]
  navigationTree.value = tree
}

function syncPathAndBuild() {
  currentPath.value = joinRoutePathParam(route.params.path)
  buildTree()
}

watch(() => route.params.path, syncPathAndBuild, { immediate: true })
watch(categoryData, () => buildTree())

onMounted(() => {
  buildTree()
})
</script>

<style scoped>
.navigation-tree {
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

.navigation-tree :deep(.tree-group + .tree-group) {
  border-top: 1px solid var(--line);
  padding-top: var(--sp-3);
}

.navigation-tree :deep(.tree-label) {
  font-size: 13px;
  font-weight: 700;
  letter-spacing: -0.01em;
  color: var(--fg);
  padding: 4px 0 var(--sp-2);
  line-height: 24px;
}

.navigation-tree :deep(.article-list) {
  list-style: none;
  padding-left: 0;
  margin: 0;
}

.navigation-tree :deep(.article-item) {
  margin: 0;
}

/* ---------- 树节点（VitePress/D3 风格：纯文字链接 + hover 变色 + 激活指示条） ---------- */
.navigation-tree :deep(.tree-item) {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  width: 100%;
  padding: 4px 0;
  font-size: 13px;
  font-weight: 500;
  line-height: 24px;
  color: var(--fg-2);
  text-align: left;
  text-decoration: none;
  border: none;
  background: transparent;
  cursor: pointer;
  position: relative;
  transition: color var(--dur-fast) ease;
}

.navigation-tree :deep(.tree-item__text) {
  flex-grow: 1;
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.navigation-tree :deep(.tree-item:hover) {
  color: var(--primary);
}

/* ---------- 层级 1：分类直属文章（无左侧导轨） ---------- */
.navigation-tree :deep(.tree-item--l1) {
  padding-left: 2px;
}

/* ---------- 目录（可折叠分支）：文字 + 右侧旋转 caret ---------- */
.navigation-tree :deep(.tree-item--folder) {
  color: var(--fg-2);
}

.navigation-tree :deep(.tree-item--folder:hover) {
  color: var(--fg);
}

.navigation-tree :deep(.tree-item--folder .tree-caret) {
  flex-shrink: 0;
  font-size: 11px;
  color: var(--fg-3);
  transition:
    color var(--dur-fast) ease,
    transform var(--dur-base) var(--ease-out);
}

.navigation-tree :deep(.tree-item--folder:hover .tree-caret) {
  color: var(--fg-2);
}

.navigation-tree :deep(.tree-caret--open) {
  transform: rotate(90deg);
}

.navigation-tree :deep(.tree-item--branch-open) {
  color: var(--fg);
}

/* ---------- 子列表导轨：左侧发丝线 + 缩进 ---------- */
.navigation-tree :deep(.tree-sublist) {
  border-left: 1px solid var(--line);
  padding-left: 16px;
}

/* ---------- 激活项：主色 + 左侧 2px 指示条（压住导轨线），不加粗 ---------- */
.navigation-tree :deep(.tree-item--active) {
  color: var(--primary);
}

.navigation-tree :deep(.tree-item--active:hover) {
  color: var(--primary);
}

.navigation-tree :deep(.tree-item--l2.tree-item--active)::before,
.navigation-tree :deep(.tree-item--l3.tree-item--active)::before {
  content: '';
  position: absolute;
  top: 6px;
  bottom: 6px;
  left: -17px;
  width: 2px;
  border-radius: 2px;
  background: var(--primary);
}
</style>
