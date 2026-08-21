<template>
  <article
    class="post-item"
    v-reveal
    :style="{ '--reveal-delay': Math.min(index, 5) * 40 + 'ms' }"
  >
    <span class="post-item__index num" aria-hidden="true">{{
      String(ordinal).padStart(2, '0')
    }}</span>

    <router-link :to="articlePath" class="post-item__main">
      <div class="post-item__meta">
        <span class="post-item__cat">{{ categoryLabel }}</span>
        <span class="divider-v"></span>
        <span class="post-item__date"><i class="fas fa-calendar-alt"></i>{{ post.date }}</span>
        <span class="divider-v"></span>
        <span class="post-item__words"
          ><i class="fas fa-file-lines"></i>{{ post.wordCount }} {{ t('wordUnit') }}</span
        >
        <span class="divider-v"></span>
        <span class="post-item__reading"
          ><i class="fas fa-clock"></i>{{ readingTime(post.wordCount) }}</span
        >
      </div>
      <h3 class="post-item__title">{{ post.title }}</h3>
      <p class="post-item__preview">{{ post.preview }}</p>
    </router-link>

    <div class="post-item__tags">
      <span
        v-for="tag in post.tags"
        :key="tag"
        class="post-item__tag"
        @click="emit('tagClick', tag)"
        >{{ tag }}</span
      >
    </div>

    <router-link :to="articlePath" class="post-item__arrow" :aria-label="post.title">
      <i class="fas fa-arrow-right"></i>
    </router-link>
  </article>
</template>

<script setup lang="ts">
import { useI18n } from 'vue-i18n'
import { readingTimeMinutes } from '@/utils/readingTime'
import type { Post } from '@/types'

interface PostCardProps {
  post: Post
  /** 列表内序号（用于 reveal 延迟） */
  index: number
  /** 全局序号（跨页连续） */
  ordinal: number
  /** 格式化后的分类标签 */
  categoryLabel: string
  /** 文章路由路径（已带 locale 前缀） */
  articlePath: string
}

defineProps<PostCardProps>()

const emit = defineEmits<{ tagClick: [tag: string] }>()

const { t } = useI18n()

function readingTime(wordCount: number) {
  return t('postReadingTime', { minutes: readingTimeMinutes(wordCount) })
}
</script>

<style scoped>
.post-item {
  display: grid;
  grid-template-columns: 64px 1fr auto;
  grid-template-areas:
    'index main arrow'
    'index tags arrow';
  align-items: center;
  gap: var(--sp-2) var(--sp-4);
  padding: var(--sp-6) 0;
  border-bottom: 1px solid var(--line);
  position: relative;
}

.post-item:last-child {
  border-bottom: none;
}

.post-item:hover {
  background: transparent;
}

.post-item__index {
  grid-area: index;
  align-self: start;
  padding-top: 6px;
  font-size: 22px;
  font-weight: 700;
  color: var(--line-strong);
  letter-spacing: 0.02em;
  font-variant-numeric: tabular-nums;
  transition: color 0.14s ease;
}

.post-item__main {
  grid-area: main;
  min-width: 0;
  display: block;
  text-decoration: none;
  color: inherit;
}

.post-item__meta {
  display: flex;
  align-items: center;
  gap: var(--sp-3);
  flex-wrap: wrap;
  font-size: var(--text-xs);
  color: var(--fg-3);
}

.post-item__meta i {
  font-size: 12px;
  color: var(--primary-muted);
  margin-right: 6px;
  transition: color 0.14s ease;
}

.post-item:hover .post-item__meta i {
  color: var(--primary);
}

.post-item__cat {
  font-weight: 600;
  color: var(--fg-2);
}

.post-item__date {
  letter-spacing: 0.01em;
  font-size: var(--text-xs);
}

.post-item__words {
  color: var(--fg-3);
}

.post-item__title {
  font-size: var(--text-xl);
  font-weight: 600;
  letter-spacing: -0.015em;
  line-height: 1.35;
  margin: var(--sp-3) 0 var(--sp-2);
  color: var(--fg);
  transition: color 0.14s ease;
}

.post-item:hover .post-item__title {
  color: var(--primary);
}

.post-item__preview {
  font-size: var(--text-md);
  color: var(--fg-2);
  line-height: 1.65;
  margin: 0;
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
  max-width: 640px;
}

.post-item__tags {
  grid-area: tags;
  display: flex;
  gap: var(--sp-2);
  margin-top: var(--sp-3);
  flex-wrap: wrap;
}

.post-item__tag {
  font-size: var(--text-xs);
  font-weight: 500;
  padding: 3px 11px;
  border-radius: var(--radius-pill);
  background: var(--surface-2);
  color: var(--fg-2);
  cursor: pointer;
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
}

.post-item__tag:hover {
  color: var(--primary);
  background: var(--tint);
}

.post-item__arrow {
  grid-area: arrow;
  align-self: center;
  display: flex;
  align-items: center;
  justify-content: center;
  color: var(--fg-3);
  font-size: 15px;
  text-decoration: none;
  padding: 6px;
  opacity: 0;
  transition:
    opacity var(--dur-base) ease,
    color var(--dur-fast) ease,
    transform var(--dur-fast) ease;
}

.post-item:hover .post-item__arrow {
  opacity: 1;
  transform: translateX(3px);
}

.post-item__arrow:hover {
  color: var(--primary);
}

@media (max-width: 767px) {
  .post-item {
    padding: var(--sp-5) 0;
    grid-template-columns: 1fr auto;
    grid-template-areas:
      'main arrow'
      'tags arrow';
  }

  .post-item__index {
    display: none;
  }

  .post-item__title {
    font-size: 17px;
  }

  .post-item__arrow {
    opacity: 1;
  }
}
</style>
